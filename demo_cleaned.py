from sympy import symbols, groebner, Poly, expand, Number, Add, factor
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
current_x = [sym for sym in x]
print(current_x)
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
        current_x[tar_bit] = current_x[tar_bit] + current_x[ctr_bits[0]]
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

##############################################
def simplify_poly(bit, output, garbage, current_x, cost_lut):
    lines = []
    cost = 0
    cost_original = 0
    if bit<len(garbage) and garbage[bit] == '-':
        print(f"output_bit: {output_bit[bit]}")
        
        f1 = current_x[bit]
        
        current_output = [sym for sym in output]
        current_output.append(f1)
        # print(output)
        # print(current_output)
        f1 = expand(f1)
        # print('initial_poly: ', f1)

        # Compute the Grobner basis
        # GB = groebner(current_output, vars, domain='GF(2)', order='lex')
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

        GB = groebner(current_output, vars, domain='GF(2)', order='lex')
        # print(f"GB: {GB}")
        # f1 = Poly(f1)
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
        
        # with open(new_file_name, 'a') as f:
        for i in GB:
            # print(f'i: {i}, type: {type(i)}')
            terms_i = Add.make_args(i)
            # print(f'terms_i: {terms_i}, type: {type(terms_i)}')
            # Process each term to reduce the power of each variable to 1
            updated_terms_i = []
            for term_i in terms_i:
                updated_term_i = term_i
                # print(f'term_i: {term_i}, type: {type(term_i)}')
                for var_i in term_i.free_symbols:
                    # print(f'type(var): {type(var)}, var: {var}')

                    term_i = Poly(term_i)
                    updated_term_i = updated_term_i / var_i**max(term_i.degree(var_i)-1, 0)
                    # print(f'updated_term_i: {updated_term_i}')
                updated_terms_i.append(updated_term_i)

            # Reconstruct the expression with reduced variable powers
            updated_expression_i = Add(*updated_terms_i)
            i = updated_expression_i

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
            # print(f"i_new: {i}")
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
            product = new_poly*i
            new_poly += i
            new_poly += product
            # print(f"new_poly_new: {new_poly}")
            new_poly = expand(new_poly.as_expr())
            
            ##############################
            # print('new_poly_i: ',new_poly)
            temp_poly = 0
            # print(f"new_poly: {new_poly}", type(new_poly))
            
            # Separate the terms in the expression
            poly_as_add = new_poly.as_expr()
            terms_new_poly = Add.make_args(poly_as_add)
            # print(terms_new_poly)
            # terms_new_poly = new_poly.all_terms()
            # print(f'terms: {terms_new_poly}, type: {type(terms_new_poly)}')
            # Process each term to reduce the power of each variable to 1
            updated_terms_new_poly = []
            for term_new_poly in terms_new_poly:
                updated_term_new_poly = term_new_poly
                # print(f'term: {term_new_poly}, type: {type(term_new_poly)}')
                for var_new_poly in term_new_poly.free_symbols:
                    # print(f'type(var): {type(var_new_poly)}, var: {var_new_poly}')

                    term_new_poly = Poly(term_new_poly)
                    updated_term_new_poly = updated_term_new_poly / var_new_poly**max(term_new_poly.degree(var_new_poly)-1, 0)
                    # print(f'updated_term: {updated_term_new_poly}')
                updated_terms_new_poly.append(updated_term_new_poly)

            # Reconstruct the expression with reduced variable powers
            updated_expression_new_poly = Add(*updated_terms_new_poly)
            new_poly = updated_expression_new_poly
            # print(f'new_poly_updated: {new_poly}')
            poly_as_add_for_j = new_poly.as_expr()
            terms_new_poly_for_j = Add.make_args(poly_as_add_for_j)
            for j in terms_new_poly_for_j:       
                # if the constant coefficient of some variable is even, then set the coeffecient to 0, else set it to 1
                # print(f"j: {j}, type: {type(j)}, len: {len(j.args)}")
                if j.is_constant():
                    if j%2 == 0:
                        new_poly = new_poly.replace(j, 0)
                        # modified_poly.append(0)
                    else:
                        new_poly = new_poly.replace(j, 1)
                        temp_poly+=1
                elif len(j.args) >0:
                    # print(f"j.args[0]: {j.args[0]}, %2: {j.args[0] % 2}")
                    if j.args[0] % 2 == 0:
                        k = j
                        for m in k.args:
                            if m.is_constant():
                                k = k.replace(m, 0)
                                break
                        new_poly = new_poly.replace(j, k)
                        temp_poly+=k
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
                        temp_poly+=k
                        # print(f"new_poly_coef_upd: {new_poly}")
                        # print(f"modified_poly: {modified_poly}")
                else:
                    temp_poly+=j
            # if i is not 0, then add a line to the list of strings that use a toffoli gate to compute the polynomial to the output variable in .tfc format
            # new_poly = Poly(dict(modified_poly), x, y, z, w, domain='ZZ')
            # print(f"modified_poly: {modified_poly}")
            # print(f'modi: {modified_poly}, type: {type(modified_poly)}')
            new_poly = temp_poly
            # new_poly = Poly(new_poly)
            # print(f"new_poly_fi: {new_poly}", type(new_poly))
            







            
        # print('new_poly: ',new_poly)
        modified_poly = 0
        # print(f"new_poly: {new_poly}", type(new_poly))
        if isinstance(new_poly, Number):
            # print("new_poly is a Number")
            if new_poly!=0:
                line = 't1 '
                line += str(output_bit[bit])
                line += '\n'
                lines.append(line)
                cost += 1
        elif len(new_poly.atoms())==1:
            # print(f"len(new_poly.atoms()): {len(new_poly.atoms())}, new_poly.atoms(): {new_poly.atoms()}")
            # print("len(new_poly.args[0].args): ", len(new_poly.args[0].args))
            line = 't2 '
            line += str(new_poly.atoms()).strip('{}')
            line += ' '
            line += str(output_bit[bit])
            line += '\n'
            lines.append(line)
            cost += 1
        else:
            # Separate the terms in the expression
            poly_as_add = new_poly.as_expr()
            terms_new_poly = Add.make_args(poly_as_add)
            # print(terms_new_poly)
            # terms_new_poly = new_poly.all_terms()
            # print(f'terms: {terms_new_poly}, type: {type(terms_new_poly)}')
            # Process each term to reduce the power of each variable to 1
            updated_terms_new_poly = []
            for term_new_poly in terms_new_poly:
                updated_term_new_poly = term_new_poly
                # print(f'term: {term_new_poly}, type: {type(term_new_poly)}')
                for var_new_poly in term_new_poly.free_symbols:
                    # print(f'type(var): {type(var_new_poly)}, var: {var_new_poly}')

                    term_new_poly = Poly(term_new_poly)
                    updated_term_new_poly = updated_term_new_poly / var_new_poly**max(term_new_poly.degree(var_new_poly)-1, 0)
                    # print(f'updated_term: {updated_term_new_poly}')
                updated_terms_new_poly.append(updated_term_new_poly)

            # Reconstruct the expression with reduced variable powers
            updated_expression_new_poly = Add(*updated_terms_new_poly)
            new_poly = updated_expression_new_poly
            # print(f'new_poly_updated: {new_poly}')
            poly_as_add_for_j = new_poly.as_expr()
            terms_new_poly_for_j = Add.make_args(poly_as_add_for_j)
            for j in terms_new_poly_for_j:       
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
                        # print(f"new_poly_coef_upd: {new_poly}")
                        # print(f"modified_poly: {modified_poly}")
                else:
                    modified_poly+=j
            # if i is not 0, then add a line to the list of strings that use a toffoli gate to compute the polynomial to the output variable in .tfc format
            # new_poly = Poly(dict(modified_poly), x, y, z, w, domain='ZZ')
            print(f"modified_poly: {modified_poly}")
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
                        lines.append(line)
                        cost += 1
                    ############################## here's the factorization part
                    # group every two terms together
                    terms_group = []
                    poly_expression = new_poly.as_expr()
                    for i in range(0, len(poly_expression.args), 2):
                        if i == len(poly_expression.args)-1:
                            terms_group.append(poly_expression.args[i])
                        else:
                            terms_group.append(poly_expression.args[i]+poly_expression.args[i+1])
                    print(terms_group)
                    for term in terms_group:

                        if term.is_Add:
                            factorized_term = factor(term)
                            print(factorized_term)
                            print(type(factorized_term))
                            can_fac = False
                            if factorized_term.is_Mul:
                                print(factorized_term.args[-1])
                                if factorized_term.args[-1].is_Add:
                                    # find a variable not in the factorized_term
                                    new_var = ''
                                    for i in range(len(output_bit)):
                                        if (output_bit[i] not in factorized_term.free_symbols) and (i != bit):
                                            new_var = output_bit[i]
                                            break
                                    print("new_var: ", new_var)

                                    if new_var != '':
                                        can_fac = True
                                        for s in range(2):
                                            for j in factorized_term.args[-1].args:
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
                                                line += new_var
                                                line += ' \n'
                                                lines.append(line)
                                                if count<=19:
                                                    cost += cost_lut[count]
                                                else:
                                                    cost += 2**count-3

                                            line = 't'
                                            line += str(len(factorized_term.args)+1)
                                            line += ' '
                                            count = 0
                                            for var in factorized_term.args[:-1]:
                                                count += 1
                                                line += str(var)
                                                line += ' '

                                            line += new_var
                                            ### add the last term
                                            line += ' '
                                            line += output_bit[bit]
                                            line += '\n'
                                            lines.append(line)
                                            if count<=19:
                                                    cost += cost_lut[count]
                                            else:
                                                cost += 2**count-3
                            if not can_fac:
                                # check if term is only a constant
                                if term.is_constant():
                                    print("term is a constant")
                                    line = 't1 '
                                    line += str(output_bit[bit])
                                    line += '\n'
                                    lines.append(line)
                                    cost += 1
                                elif len(term.free_symbols) == 1:
                                    symbol = next(iter(term.free_symbols))
                                    print("Symbol in the polynomial term:", symbol)
                                    line = 't1 '
                                    line += str(symbol)
                                    line += ' '
                                    line += str(output_bit[bit])
                                    line += '\n'
                                    lines.append(line)
                                    cost += 1
                                else:
                                    for j in (term.as_expr()).args:
                                        print("j: ", j)
                                        count = 0
                                        line = 't'
                                        temp_line = ''
                                        print("len(j.args): ", len(j.args))
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
                                        print(f"OOOoutputBit[bit]: {output_bit[bit]}")
                                        line += output_bit[bit]
                                        line += '\n'
                                        lines.append(line)
                                        if count<=19:
                                            cost += cost_lut[count]
                                        else:
                                            cost += 2**count-3

                        else:
                            print("term_not_add: ", term)
                            term = term.as_expr()
                            count = 0
                            line = 't'
                            temp_line = ''
                            print("len(term.args): ", len(term.args))
                            if len(term.args) != 0:
                                for k in term.args:
                                    count += 1
                                    temp_line += str(k)
                                    temp_line += ' '
                                # count += 1
                            elif term != 1:
                                temp_line += str(term)
                                temp_line += ' '
                                count += 1
                            line += str(count+1)
                            line += ' '
                            line += temp_line
                            # line += ' '
                            print(f"OOOoutputBit[bit]: {output_bit[bit]}")
                            line += output_bit[bit]
                            line += '\n'
                            lines.append(line)
                            if count<=19:
                                cost += cost_lut[count]
                            else:
                                cost += 2**count-3




                    ##############################    
                    # for j in (new_poly.as_expr()).args:
                    #     count = 0
                    #     line = 't'
                    #     temp_line = ''
                    #     # print(len(j.args))
                    #     if len(j.args) != 0:
                    #         for k in j.args:
                    #             count += 1
                    #             temp_line += str(k)
                    #             temp_line += ' '
                    #         # count += 1
                    #     elif j != 1:
                    #         temp_line += str(j)
                    #         temp_line += ' '
                    #         count += 1
                    #     line += str(count+1)
                    #     line += ' '
                    #     line += temp_line
                    #     # line += ' '
                    #     # print(f"outputBit[bit]: {output_bit[bit]}")
                    #     line += output_bit[bit]
                    #     line += '\n'
                    #     lines.append(line)
                    #     if count<=19:
                    #         cost += cost_lut[count]
                    #     else:
                    #         cost += 2**count-3
                    ##############################
                            
                elif new_poly.degree()==0 and new_poly!=0:
                    line = 't1 '
                    line += str(output_bit[bit])
                    line += '\n'
                    lines.append(line)
                    cost += 1
                    # print(f"new_poly is a Poly: {new_poly}, cost: {cost}")
            elif isinstance(new_poly, Number):
                if new_poly!=0:
                    line = 't1 '
                    line += str(output_bit[bit])
                    line += '\n'
                    lines.append(line)
                    cost += 1
                    # print(f"new_poly is a Number: {new_poly}, cost: {cost}")
    # lines.append('\n')
    return lines, cost, cost_original    
    


##############################################
from joblib import Parallel, delayed
results = []
times = numbits
num_parallel = 10
for bit in range(0, times, num_parallel):
    result = Parallel(n_jobs=-1)(delayed(simplify_poly)(j, output, garbage, current_x, cost_lut) for j in range(bit, bit+num_parallel))
    results.extend(result)
    print(f"bit: {bit}")

lines, costs, costs_original = zip(*results)

# print(f"lines: {lines}")
# print(f"costs: {costs}")

with open(new_file_name, 'a') as f:
    for lines_sub in lines:
        for line in lines_sub:
            f.write(line)

for cost in costs:
    cost_tot += cost
for cost_original in costs_original:
    cost_org_tot += cost_original
print('initial circuit cost: ', cost_init)
print('original total cost: ', cost_org_tot)
print('minimized total cost: ', cost_tot)

