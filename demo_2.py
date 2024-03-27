from sympy import symbols, groebner, Poly, expand, Number, Add
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-REAL', type=str, 
                    help='Input file (.real format)')
parser.add_argument('-output', type=str, 
                    help='Output file (.real format)')
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
      var_real = []
      for idx, bit_name in enumerate(rest.split(' ')):
        bit_dict[bit_name] = idx
        var_real.append(bit_name)
    elif label == 'inputs':
        input_bit = []
        for idx, bit_name in enumerate(rest.split(' ')):
            input_bit.append(bit_name)
    elif label == 'outputs':
      output_bit = []
      for idx, bit_name in enumerate(rest.split(' ')):
        output_bit.append(bit_name)
    elif label == 'constants':
       constants = rest
    elif label == 'garbage':
       garbage_original = list(rest)
       garbage = list(reversed(rest))

gates = [line.rstrip() for line in content if line[0] != '.']
print('  - numgates:\t{}'.format(len(gates)))

new_file_name = cmd_args.output
with open(new_file_name, 'w') as f:
    f.write(f'.numvars {numbits}\n')
    f.write(f'.variables {" ".join(bit_dict.keys())}\n')
    f.write(f'.inputs {" ".join(input_bit)}\n')
    f.write(f'.outputs {" ".join(output_bit)}\n')
    f.write(f'.constants {"".join(constants)}\n')
    # f.write(f'.garbage {"".join(garbage)}\n')
    f.write(f'.garbage {"".join(garbage_original)}\n')
    # f.write(f'.numgates {len(gates)}\n')
    # for gate in gates:
    #     f.write(f'{gate}\n')
f.close()
# x = [symbols(('x'+str(i))) if constants[i] == '-' else int(constants[i]) for i in range(numbits)]
x = [symbols(var_real[i]) if constants[i] == '-' else int(constants[i]) for i in range(numbits)]
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
        print(f"output_bit: {output_bit[bit]}")
        cost = 0
        f1 = current_x[bit]
        cost_original = 0
        current_output = [sym for sym in output]
        current_output.append(f1)
        # print(output)
        # print(current_output)
        f1 = expand(f1)
        print('initial_poly: ', f1)
        # Compute the Grobner basis
        GB = groebner(current_output, vars, order='lex')
        # 'lex' (lexicographic order)
        # 'grlex' (graded lexicographic order)
        # 'grevlex' (graded reverse lexicographic order)
        # print(GB)

        # add all polynomials together to make a new polynomial
        new_poly = 0

        #calculate the cost of f1
#     # Separate the terms in the expression
        terms = Add.make_args(f1)
        # print(f'terms: {terms}, type: {type(terms)}')
        # Process each term to reduce the power of each variable to 1
        updated_terms = []
        for term in terms:
            updated_term = term
            # print(f'term: {term}, type: {type(term)}')
            for var in term.free_symbols:
                # print(f'type(var): {type(var)}, var: {var}')

                term = Poly(term)
                updated_term = updated_term / var**max(term.degree(var)-1, 0)
                # print(f'updated_term: {updated_term}')
            updated_terms.append(updated_term)

        # Reconstruct the expression with reduced variable powers
        updated_expression = Add(*updated_terms)
        f1 = updated_expression

        # print("Original Expression:", expression)
        # print("Expression with Reduced Variable Powers:", updated_expression)
        # print("f1_new: ", f1)
            # if j.is_Pow:
            #     print(f'j is a power: {j}')
            #     if j.args[1] > 1:
            #         k = j
            #         for m in k.args:
            #             if m.is_constant():
            #                 k = k.replace(m, 1)
            #                 break
            #         f1 = f1.replace(j, k)
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
        # print(f"k1: {f1}")
        # new_curr_out = []
        # new_curr_out.append(f1)
        # new_curr_out += current_output
        # print(f"new_curr_out: {new_curr_out}")
        # current_output.append(f1)
        # print(f"current_output: {current_output}")
        # GB = groebner(current_output, vars, order='lex')
        # print(f"GB: {GB}")
        f1 = Poly(f1)
        # print(type(f1))
        if isinstance(f1, Poly):
            # print("f1 is a Poly")
            if f1 != 0:
                f1 = Poly(f1)
                # print(f1)
                # print(len(f1.args[0].args))
                # for every cube in the polynomial, count the number of literals in each cube
                if len(f1.free_symbols) == 1 and len(f1.args[0].args) <= 1:
                        symbol = next(iter(f1.free_symbols))
                        # print("Symbol in the polynomial:", symbol)
                        cost_original += 1
                for j in f1.args[0].args:
                    count = 0
                    if len(j.args) != 0:
                        for k in j.args:
                            count += 1
                    count += 1
                    # print(f'count: {count}')
                    if count<=19:
                        cost_original += cost_lut[count]
                        # print(f"cost_lut[count]: {cost_lut[count]}, cost_original: {cost_original}")
                    else:
                        cost_original += 2**count-3
        elif isinstance(f1, Number):
            if f1!=0:
                cost_original += 1
                # print(f"f1 is a Number: {f1}, cost: {cost_original}")
        # check the polynomials in the Grobner basis, if there is a variable with power greater than 1, then reduce the power to 1
        with open(new_file_name, 'a') as f:
            for i in GB:
                # print(f'i: {i}, type: {type(i)}')
                terms_i = Add.make_args(i)
                # print(f'terms: {terms}, type: {type(terms)}')
                # Process each term to reduce the power of each variable to 1
                updated_terms_i = []
                for term_i in terms_i:
                    updated_term_i = term_i
                    # print(f'term: {term}, type: {type(term)}')
                    for var_i in term_i.free_symbols:
                        # print(f'type(var): {type(var)}, var: {var}')

                        term_i = Poly(term_i)
                        updated_term_i = updated_term_i / var_i**max(term_i.degree(var_i)-1, 0)
                        # print(f'updated_term: {updated_term}')
                    updated_terms_i.append(updated_term_i)

                # Reconstruct the expression with reduced variable powers
                updated_expression_i = Add(*updated_terms_i)
                i = updated_expression_i
                # for j in i.args:
                #     if j.is_Pow:
                #         if j.args[1] > 1:
                #             k = j
                #             for m in k.args:
                #                 if m.is_constant():
                #                     k = k.replace(m, 1)
                #                     break
                #             i = i.replace(j, k)

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
            modified_poly = 0
            # print(f"new_poly: {new_poly}", type(new_poly))
            if isinstance(new_poly, Number):
                print("new_poly is a Number")
                if new_poly!=0:
                    line = 't1 '
                    line += str(output_bit[bit])
                    line += '\n'
                    f.write(line)
                    cost += 1
            elif len(new_poly.atoms())==1:
                print(f"len(new_poly.atoms()): {len(new_poly.atoms())}, new_poly.atoms(): {new_poly.atoms()}")
                # print("len(new_poly.args[0].args): ", len(new_poly.args[0].args))
                line = 't1 '

                cost += 1
            else:
                for j in new_poly.args[0].args:       
                    # if the constant coefficient of some variable is even, then set the coeffecient to 0, else set it to 1
                    # print(f"j: {j}, type: {type(j)}, len: {len(j.args)}")
                    if j.is_constant():
                        if j%2 == 0:
                            new_poly = new_poly.replace(j, 0)
                            # modified_poly.append(0)
                        else:
                            new_poly = new_poly.replace(j, 1)
                            modified_poly+=1
                    elif len(j.args) >0:
                        # print(f"j.args[0]: {j.args[0]}, %2: {j.args[0] % 2}")
                        if j.args[0] % 2 == 0:
                            k = j
                            for m in k.args:
                                if m.is_constant():
                                    k = k.replace(m, 0)
                                    break
                            new_poly = new_poly.replace(j, k)
                            modified_poly+=k
                        else:
                            k = j
                            for m in k.args:
                                if m.is_constant():
                                    # print(f"m: {m}")
                                    k = k.replace(m, 1)
                                    # print(f"k: {k}")
                                    break
                            new_poly = Poly(new_poly)
                            new_poly = new_poly.replace(j, k)
                            modified_poly+=k
                            # print(f"new_poly: {new_poly}")
                    else:
                        modified_poly+=j
                # if i is not 0, then add a line to the list of strings that use a toffoli gate to compute the polynomial to the output variable in .tfc format
                # new_poly = Poly(dict(modified_poly), x, y, z, w, domain='ZZ')
                # print(f"modified_poly: {modified_poly}")
                # print(f'modi: {modified_poly}, type: {type(modified_poly)}')
                new_poly = modified_poly
                new_poly = Poly(new_poly)
                if isinstance(new_poly, Poly):
                    # print("new_poly is a Poly")
                    # print("new_poly.degree(): ", new_poly.degree())
                    if not new_poly.degree() == 0:
                        # new_poly = Poly(new_poly)
                        # print(new_poly)
                        # for every cube in the polynomial, count the number of literals in each cube
                        # print(new_poly.args[0])
                        # print(len(new_poly.args[0].args))
                        if len(new_poly.free_symbols) == 1 and len(new_poly.args[0].args) <= 1:
                            symbol = next(iter(new_poly.free_symbols))
                            # print("Symbol in the polynomial:", symbol)
                            line = 't1 '
                            line += str(symbol)
                            line += ' '
                            line += str(output_bit[bit])
                            line += '\n'
                            f.write(line)
                            cost += 1
                        for j in new_poly.args[0].args:
                            count = 0
                            line = 't'
                            temp_line = ''
                            # print(len(j.args))
                            if len(j.args) != 0:
                                for k in j.args:
                                    count += 1
                                    temp_line += str(k)
                                    temp_line += ' '
                                # count += 1
                            elif j != 1:
                                temp_line += str(j)
                                temp_line += ' '
                                count += 1
                            line += str(count+1)
                            line += ' '
                            line += temp_line
                            # line += ' '
                            # print(f"outputBit[bit]: {output_bit[bit]}")
                            line += output_bit[bit]
                            line += '\n'
                            f.write(line)
                            if count<=19:
                                cost += cost_lut[count]
                            else:
                                cost += 2**count-3
                    elif new_poly.degree()==0 and new_poly!=0:
                        line = 't1 '
                        line += str(output_bit[bit])
                        line += '\n'
                        f.write(line)
                        cost += 1
                        # print(f"new_poly is a Poly: {new_poly}, cost: {cost}")
                elif isinstance(new_poly, Number):
                    if new_poly!=0:
                        line = 't1 '
                        line += str(output_bit[bit])
                        line += '\n'
                        f.write(line)
                        cost += 1
                        # print(f"new_poly is a Number: {new_poly}, cost: {cost}")
        f.close()
        print(new_poly)
        print(f"cost: {cost}")
        print(f"cost_original: {cost_original}")        
        cost_tot += cost
        cost_org_tot += cost_original

print('initial circuit cost: ', cost_init)
print('original total cost: ', cost_org_tot)
print('minimized total cost: ', cost_tot)



