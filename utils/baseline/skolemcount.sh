#!/usr/bin/bash
# set -x

python="python"

function stop()
{
  echo "c it $((it))"
  echo "s fc 2 ** $skolemcountest"
  echo "c time $SECONDS sec"
  exit 0
}
trap 'stop' SIGINT


usage () {
cat <<EOF
usage: $0 <options> <filename.qdimacs>
where <option> is one of the following:
  -h, --help        print this message and exit
EOF
  exit 0
}

die () {
  echo "*** skolemcount.sh: $*" 1>&2
  exit 1
}

msg () {
  echo "[skolemcount.sh] $*"
}

set_env () {
  filename=$opt
  samplefile=samples.out
  update_nos=0
  filecount=0
  unisampseed=345
  samplenumest=0


  basefile=$(basename $filename)
  if [ $basefile != $filename ] ; then
    cp $filename ./
    filename=$basefile
  fi
  if [ $absolute_count == false ] ; then
    echo "[skolemfc] counting $filename with epsilon = $epsilon, delta = $delta"
  else
    echo "[skolemfc] counting $filename with Baseline"
  fi
  foldername=${filename%%.qdimacs}-files
  cnffile=${filename%%.qdimacs}.cnf

  num_vars=$(grep "^p" $filename | awk '{print  $3}')
  num_clauses=$(grep "^p" $filename | awk '{print  $4}')
  num_e_vars=`expr $(grep "^e" $filename | awk '{ print NF }') - 2`
  num_a_vars=`expr $(grep "^a" $filename | awk '{ print NF }') - 2`

  counterrestext="^c s exact"
  counterrespos=6
  num_new_clauses=`expr $num_clauses + $num_vars`
  samplingcnffile="sample_F2_${filename}.cnf"
  samplingxcnffile="sample_F2_${filename}_gpmc.cnf"

  ecounter="./bin/gpmc"

  echo "[skolemfc] vars = $num_vars (a $num_a_vars + e $num_e_vars) clauses = $num_clauses"
}

create_sample_or_allsol () {
  ## Create the necessary files
  samplingcnffile="sample_F2_${filename}.cnf"
  aline=`(grep "^a" $filename)`
  echo "${aline//a/c ind}" > $cnffile
  grep -P '^(?!(a|e))' $filename >> $cnffile

  $python f2_formula_creator.py $filename -x

  if [ $get_baseline_samps == true ] ; then
    echo "[skolemfc] counting number of solutions for sampling formula"
    num_sol_f2=$($pcounter $samplingcnffile > >(tee counting.log) 2>&1 | grep "^s mc" | awk '{print  $3}')
    echo "[skolemfc] Counter calls needed for Baseline $num_sol_f2"
    exit 0
  fi

  if [[ $samplenumest -eq 0 ]] && [[ $absolute_count == false ]] ;then
      get_est_sample_stoppingrule;
  fi

  num_counts_needed=$samplenumest

  if [ $absolute_count == true ] ; then
    echo "[skolemfc] employing cryptominisat to generate all solutions"
    ./bin/cryptominisat5 --maxsol 2 --verb 0 --dumpresult $samplefile --onlysampling $samplingcnffile > sampling.log
    num_counts_needed=`cat $samplefile | wc -l`
    echo "[skolemfc] cryptominisat generated all ($num_counts_needed) solutions of the formula ($SECONDS s)"
    num_sol_f2=$num_counts_needed
  else
    echo "[skolemfc] employing unisamp to generate $num_counts_needed samples [num_sol : $num_sol_f2]"
    timeout $timeout ./bin/unisamp --seed $unisampseed --samples $num_counts_needed --epsilon $epsilons --sampleout $samplefile  $samplingcnffile > sampling.log
    if [ $timeout -le $SECONDS ]; then
      echo "[skolemfc] timeout reached, exiting from unisamp calls"
      stop
    fi
    # echo -e "[skolemfc] unisamp generated $num_counts_needed samples from $samplingcnffile \e[0;33m  ($SECONDS s) \e[0m"
    echo "[skolemfc] unisamp generated $num_counts_needed samples from $samplingcnffile ($SECONDS s)"
  fi
}

gen_files () {
  mkdir $foldername
  while read -r line
  do
    ithfile=$foldername/$filecount-$cnffile
    cp $cnffile $ithfile
    echo "p cnf $num_vars $num_new_clauses" > $ithfile
    grep -P '^(?!(p|c))' $cnffile >> $ithfile
    for word in $line; do
      if [ $word != "0" ]; then
        echo "$word 0" >> $ithfile
      fi
    done
    filecount=`expr $filecount + 1`
    if [ $timeout -le $SECONDS ]; then
      echo "[skolemfc] timeout reached, exiting in gen_files after $filecount files"
      stop
    fi
  done < <(cat $samplefile)
  echo "[skolemfc] generated $filecount files for counting ($SECONDS s)"
}

count_absolute () {
  it=1
  (( filenum = it - 1 ))
  sum_count=0
  echo "[skolemfc] Counting Starts Here" > counting.log
  echo "[skolemfc] Iteration |         count now          | perc reached | time elapsed"
  while [[ $it -le $filecount ]] && [[ $SECONDS -le $timeout ]] ; do
    ithfile=$foldername/$filenum-$cnffile
    echo "[skolemfc] Counting file ${it}" >> counting.log
    timeremaining=$((timeout - SECONDS))
    # check if ${ithfile} exists
    ls ${ithfile} > /dev/null 2>&1
    if [ $? -ne 0 ]; then
      echo "[skolemfc] file ${ithfile} does not exist, serious error"
      stop;
    fi
    expcount=$(timeout $timeremaining $ecounter $ithfile > >(tee counting.log) 2>&1 | grep "^c s exact arb int " | awk '{print  $6}')
    count=`echo "scale=5; l ($expcount) / l (2)" | bc -l`

    sum_count=`echo "( $sum_count + $count ) " | bc`
    skolemcountest=`echo "scale=2; ( $sum_count * $filecount / $it ) " | bc`
    if (( $it % 1 == 0 )); then
      conf_now=`echo  "scale=4; ( $it  * 100 /  $filecount  )" | bc -l `
      echo "[skolemfc]      $it |  $skolemcountest  |   $conf_now    | $SECONDS s"
    fi
    ((it = it + 1))
  done
  ((it = it - 1))
  skolemcountest=$sum_count
}


absolute_count=true
limit=0
get_baseline_samps=false
exactmode=true
timeout=36000

while [ $# -gt 0 ]
do
  opt=$1
  case $opt in
    -h|--help) usage;;
  esac
  shift
done


set_env;
create_sample_or_allsol;
gen_files;
count_absolute;

confindent="conf"
stop;

