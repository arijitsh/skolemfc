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

get_est_sample_stoppingrule () {
  # create sigma
  ./bin/unisamp --seed 345 --samples 1 --sampleout $samplefile  $samplingcnffile > sampling.log
  # create F ^ (X = \sigma) and count it
  gen_files;
  ithfile=$foldername/1-$cnffile
  expcount=$(./bin/appmcplus4s2 $ithfile > >(tee counting.log) 2>&1 | grep "^s mc " | awk '{print  $3}')
  two="2"
  samplenumest=`echo "$stoppingruletarget * l ($two) / (l ($expcount)) " | bc -l` #TODO remove 1.2 if possible
  #ceiling function
  samplenumest=${samplenumest%%.*}
  echo "[skolemfc] Estimate Number of samples = $samplenumest ($num_e_vars / log ($expcount) )"
}

set_env () {
  filename=$opt
  samplefile=samples.out
  update_nos=0
  filecount=0
  unisampseed=345
  samplenumest=0

  epsilonf=`echo "scale=4; ( $epsilon * 0.6 ) " | bc -l `
  deltaf=`echo "scale=4; ( $delta * 0.5 ) " | bc -l `
  epsilons=`echo "scale=4; ( $epsilon * 0.3 ) " | bc -l `

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

  stoppingruletarget=`echo  "(4.0 * $num_e_vars *  l (2 / $deltaf) * (1 + $epsilonf) /( $epsilonf * $epsilonf) )" | bc -l `

  num_counts_needed=`echo  "(100.0 * $num_e_vars *  l (2 / $delta) + 1)" | bc -l `
  num_counts_needed=${num_counts_needed%%.*} ## ceiling function

  counterrestext="^c s exact"
  counterrespos=6
  pmode="-x"
  num_new_clauses=`expr $num_clauses + $num_vars`
  appmcdelta=`echo "1 / $num_counts_needed" | bc -l`
  samplingcnffile="sample_F2_${filename}.cnf"
  samplingxcnffile="sample_F2_${filename}_gpmc.cnf"

  counter="./bin/appmcplus4s2 --pivotbysqrt2 1 --d $appmcdelta"
  pcounter="./bin/gpmc -mode=2"
  ecounter="./bin/ganak"

  echo "[skolemfc] vars = $num_vars (a $num_a_vars + e $num_e_vars) clauses = $num_clauses"
  $python f2_formula_creator.py $filename -a

  if [ $num_counts_needed -ge $limit ] && [ $limit -ne 0 ]  ; then
    actual_counts_needed=$num_counts_needed
    num_counts_needed=$limit
    echo "[skolemfc] needed $actual_counts_needed many counts, capping to $num_counts_needed counts for simplicity"
  fi

  if [ $exactmode == false ] ; then
    setapproximation;
  fi

  if [ $absolute_count == false ] ; then
    $python f2_formula_creator.py $filename $pmode
    echo "[skolemfc] GPMC counting the number of solutions for sampling formula <$samplingcnffile>"
    num_sol_f2=$($pcounter $samplingxcnffile > >(tee counting.log) 2>&1 | grep "$counterrestext" | awk -v field=$counterrespos '{print $field}' )
    if [ $num_sol_f2 == 0 ] ; then
      echo "G formula has no solutions"
      echo "s amc 0 conf"
      exit 0
    fi
    echo "[skolemfc] GPMC returned the count $num_sol_f2"
  fi
}

create_sample_or_allsol () {
  ## Create the necessary files
  samplingcnffile="sample_F2_${filename}.cnf"
  aline=`(grep "^a" $filename)`
  echo "${aline//a/c ind}" > $cnffile
  grep -P '^(?!(a|e))' $filename >> $cnffile

  $python f2_formula_creator.py $filename $pmode

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
    ./bin/cryptominisat5_allsol --maxsol 2 --verb 0 --dumpresult $samplefile --onlysampling $samplingcnffile > sampling.log
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

setapproximation(){
  rm -rf $foldername
  ecounter=$counter
  filecount=0
  pcounter=$counter
  counterrestext="^s mc"
  counterrespos=3
  pmode="-a"
  samplingxcnffile="sample_F2_${filename}.cnf"
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
    expcount=$(timeout $timeremaining $ecounter $ithfile > >(tee counting.log) 2>&1 | grep "^s mc " | awk '{print  $3}')
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


epsilon=0.8
delta=0.8
absolute_count=true
only_count=false
only_sample=false
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

