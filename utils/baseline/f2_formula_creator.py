import sys

print("[skolemfc->formula_creator] got formula %s to create sampling formula"%(sys.argv[1]))
# 1. Read all variables

constraints = []

with open(sys.argv[1],'r') as f:
    for line in f:
        words = line.split()
        if len(words) == 0:
            pass
        elif str.isdigit(words[0]) or str.isdigit(words[0][1:]):
            constraint = [eval(i) for i in words]
            assert(constraint[-1:] == [0])
            constraints.append(constraint[:-1])
        elif words[0] == 'a':
            a_vars = [eval(i) for i in words[1:-1]]
            set_a_vars = set(a_vars)
        elif words[0] == 'e':
            e_vars = [eval(i) for i in words[1:-1]]
            set_e_vars = set(e_vars)
        elif words[0] == 'p' and words[1] == "cnf":
            dec_num_vars = eval(words[2])
            dec_num_cls = eval(words[3])

print("[skolemfc->formula_creator] input formula has %d vars and %d clauses"%(dec_num_vars, dec_num_cls))

# 2. Create variable mapping

mapping = [0]
offset = dec_num_vars
missing_vars = []

for i in range(dec_num_vars):
    if (i+1) in set_a_vars:
        mapping.append(i+1)
    elif (i+1) in set_e_vars:
        offset += 1
        mapping.append(offset)
    else:
        missing_vars.append(i+1)
        mapping.append(i+1)

if missing_vars:
    print("[skolemfc] These variables are neither in e or a : ", end = '')
    print(missing_vars)

# 3. Generate New Constraits

new_consts = []

# 3a. F(X,Y')

for constraint in constraints:
    const_vars = set([abs(i) for i in constraint])
    set_const_vars = set(const_vars)
    if not set_const_vars.issubset(set_a_vars):
        new_const = []
        for lit in constraint:
            new_var = mapping[abs(lit)]
            if lit < 0:
                new_var *= -1
            new_const.append(new_var)
        new_consts.append(new_const)

# 3b. (Y != Y')
k_const = []
for var in range(dec_num_vars):
    if not (var+1 == mapping[var+1]):
        offset += 1
        const = [var+1, mapping[var+1], -1*offset]
        new_consts.append(const)
        const = [-1*(var+1), -1*mapping[var+1], -1*offset]
        new_consts.append(const)
        k_const.append(offset)
new_consts.append(k_const)

# 4. Print the new file
outfile = "sample_F2_" + sys.argv[1]
if sys.argv[2] == '-a':
    outfile = outfile + ".cnf"
elif sys.argv[2] == '-x':
    outfile = outfile + "_gpmc.cnf"
else:
    print("[skolemfc->formula_creator] Error: Expected argument \
        -a or -x after filename. got %s. EXITING"%(sys.argv[2]))
    exit(1)
num_new_const = len(constraints) + len(new_consts)
num_new_vars = offset
constraints += new_consts

with open(outfile ,'w') as f:
    f.write("p cnf " + str(num_new_vars) + " " + str(num_new_const) +'\n')
    if sys.argv[2] == '-a':
        f.write("c ind ")
    elif sys.argv[2] == '-x':
        f.write("c p show ")
    for var in a_vars:
        f.write(str(var) + " ")
    f.write("0\n")
    for const in constraints:
        for var in const:
            f.write(str(var) + " ")
        f.write("0\n")
print("[skolemfc->formula_creator] created sampling formula %s with %d vars and %d clauses"%(outfile, num_new_vars,num_new_const))
