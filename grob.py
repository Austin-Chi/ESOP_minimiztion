from sympy import symbols, groebner, Poly, expand, Number
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('REAL', type=str, 
                    help='Input file (.real format)')
cmd_args = parser.parse_args()
real_file = cmd_args.REAL
circ_name = os.path.basename(real_file).split('.real')[0]
qasm_file = circ_name + '.qasm'

realfile = open(real_file,'r', encoding="utf8")

content = realfile.readlines()
content = [line.replace('\n', '') for line in content]
content = [line for line in content if not line == '']
content = [line for line in content if line[0] != '#']

constants = []
garbage = []
cost_init = 0
cost_lut0 = [1, 1, 5, 13, 29, 61, 125, 253, 509, 1021]
for i in range(11, 20):
    cost_lut0.append(2**i-3)

for info in [line for line in content if line[0] == '.'][:-2]:
    infos = [i for i in info.split(' ') if len(i) > 0]
    label = infos[0][1:]
    rest = ' '.join(infos[1:])
    print('  - {}:\t{}'.format(label, rest))

    if label == 'numvars':
      numbits = int(rest)
      constants = ['-' for i in range(numbits)]
      garbage = ['-' for i in range(numbits)]
    elif label == 'variables':
      bit_dict = {}
      for idx, bit_name in enumerate(rest.split(' ')):
        bit_dict[bit_name] = idx
    elif label == 'constants':
       constants = rest
    elif label == 'garbage':
       garbage = list(reversed(rest))

gates = [line.rstrip() for line in content if line[0] != '.']
print('  - numgates:\t{}'.format(len(gates)))


x = [symbols(('x'+str(i))) if constants[i] == '-' else int(constants[i]) for i in range(numbits)]
vars = []
print(constants)
for i in range(numbits):
    if constants[i] == '-':
        vars.append(x[i])
print(vars)
current_x = x
output = []
for i in range(numbits):
    if constants[i] == '-':
        output.append(x[i]**2+x[i])
# print(output)

print('\nConverting...', end=' ')
for idx, g in enumerate(gates):
    g_split = g.split(' ')
    numctrl = int(g_split[0][1:]) - 1
    ctr_labels = g_split[1:-1]
    tar_label = g_split[-1]
    ctr_bits = [bit_dict[c] for c in ctr_labels]
    tar_bit = bit_dict[tar_label]
    cost_init += cost_lut0[numctrl]
    if numctrl == 0:
        current_x[tar_bit] = current_x[tar_bit] + 1
    elif numctrl == 1:
        current_x[tar_bit] = current_x[tar_bit] + x[ctr_bits[0]]
    else:
        product = 1
        for bit in ctr_bits:
           product = product * current_x[bit]
        current_x[tar_bit] = current_x[tar_bit] + product
print('circuit done.')
# print(current_x)
# Define the variables using an array of strings
# variables = ['x', 'y', 'z']
# Define the variables using the sympy librar


# the quantum cost
cost_lut = [0, 1, 1, 5, 13, 29, 61, 125, 253, 509, 1021]
for i in range(11, 20):
    cost_lut.append(2**i-3)
cost_tot = 0
cost_org_tot = 0


# Define the polynomials
for bit in range(numbits):
    if garbage[bit] == '-':
        cost = 0
        f1 = current_x[bit]
        cost_original = 0
        current_output = [sym for sym in output]
        current_output.append(f1)
        print(output)
        print(current_output)
        f1 = expand(f1)
        print('initial_poly: ', f1)
        # Compute the Grobner basis
        GB = groebner(current_output, vars, order='lex')
        # 'lex' (lexicographic order)
        # 'grlex' (graded lexicographic order)
        # 'grevlex' (graded reverse lexicographic order)
        print(GB)

        # add all polynomials together to make a new polynomial
        new_poly = 0

        #calculate the cost of f1
        for j in f1.args:
            if j.is_Pow:
                if j.args[1] > 1:
                    k = j
                    for m in k.args:
                        if m.is_constant():
                            k = k.replace(m, 1)
                            break
                    f1 = f1.replace(j, k)
        for j in f1.args:
            if j.is_constant():
                if j%2 == 0:
                    f1 = f1.replace(j, 0)
                else:
                    f1 = f1.replace(j, 1)
            elif len(j.args) >0:
                if j.args[0] % 2 == 0:
                    k = j
                    for m in k.args:
                        if m.is_constant():
                            k = k.replace(m, 0)
                            break
                    f1 = f1.replace(j, k)
                else:
                    k = j
                    for m in k.args:
                        if m.is_constant():
                            k = k.replace(m, 1)
                            break
                    f1 = f1.replace(j, k)
        if f1 != 0:
            f1 = Poly(f1)
            # for every cube in the polynomial, count the number of literals in each cube
            for j in f1.args[0].args:
                count = 0
                if len(j.args) != 0:
                    for k in j.args:
                        count += 1
                count += 1
                if count<=10:
                    cost_original += cost_lut[count]
                else:
                    cost_original += 12*count-34
        # check the polynomials in the Grobner basis, if there is a variable with power greater than 1, then reduce the power to 1
        for i in GB:
            for j in i.args:
                if j.is_Pow:
                    if j.args[1] > 1:
                        k = j
                        for m in k.args:
                            if m.is_constant():
                                k = k.replace(m, 1)
                                break
                        i = i.replace(j, k)

            for j in i.args:       
                # if the constant coefficient of some variable is even, then set the coeffecient to 0, else set it to 1
                if j.is_constant():
                    if j%2 == 0:
                        i = i.replace(j, 0)
                    else:
                        i = i.replace(j, 1)
                elif len(j.args) >0:
                    if j.args[0] % 2 == 0:
                        k = j
                        for m in k.args:
                            if m.is_constant():
                                k = k.replace(m, 0)
                                break
                        i = i.replace(j, k)
                    else:
                        k = j
                        for m in k.args:
                            if m.is_constant():
                                k = k.replace(m, 1)
                                break
                        i = i.replace(j, k)
            # if i is not 0, then add a line to the list of strings that use a toffoli gate to compute the polynomial to the output variable in .tfc format
            if not i.is_constant():
                i = Poly(i)
                # for every cube in the polynomial, count the number of literals in each cube
                for j in i.args[0].args:
                    count = 0

                    if len(j.args) != 0:
                        for k in j.args:
                            count += 1
                    count += 1
                    # if count<=10:
                    #     cost += cost_lut[count]
                    # else:
                    #     cost += 12*count-34
            # elif i.is_constant() and i!=0:
            #     cost += 1
            new_poly += i

        for j in new_poly.args:       
            # if the constant coefficient of some variable is even, then set the coeffecient to 0, else set it to 1
            if j.is_constant():
                if j%2 == 0:
                    new_poly = new_poly.replace(j, 0)
                else:
                    new_poly = new_poly.replace(j, 1)
            elif len(j.args) >0:
                if j.args[0] % 2 == 0:
                    k = j
                    for m in k.args:
                        if m.is_constant():
                            k = k.replace(m, 0)
                            break
                    new_poly = new_poly.replace(j, k)
                else:
                    k = j
                    for m in k.args:
                        if m.is_constant():
                            k = k.replace(m, 1)
                            break
                    new_poly = new_poly.replace(j, k)
        # if i is not 0, then add a line to the list of strings that use a toffoli gate to compute the polynomial to the output variable in .tfc format
        if isinstance(new_poly, Poly):
            if not new_poly.degree() == 0:
                new_poly = Poly(new_poly)
                # for every cube in the polynomial, count the number of literals in each cube
                for j in new_poly.args[0].args:
                    count = 0

                    if len(j.args) != 0:
                        for k in j.args:
                            count += 1
                    count += 1
                    if count<=10:
                        cost += cost_lut[count]
                    else:
                        cost += 12*count-34
            elif new_poly.degree()==0 and new_poly!=0:
                cost += 1
        elif isinstance(new_poly, Number):
            if new_poly !=0:
                cost +=1
        print(new_poly)
        print('current_output: ', bit)
        print(f"cost: {cost}")
        print(f"cost_original: {cost_original}")
        cost_tot += cost
        cost_org_tot += cost_original

print('initial circuit cost: ', cost_init)
print('original total cost: ', cost_org_tot)
print('minimized total cost: ', cost_tot)



