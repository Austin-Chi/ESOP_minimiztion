from sympy import symbols, groebner, Poly, expand, Number, Add

# Define the variables using an array of strings
# variables = ['x', 'y', 'z']
# Define the variables using the sympy library
x, y, z, w = symbols('x y z w')

# define output varible name
output = 'f'


# the quantum cost
cost_lut = [0, 1, 1, 5, 13, 26, 38, 50, 62, 74, 86]
cost = 0



# Define the polynomials
f1 =x + y*z*(x*w + y*z*(x + w)) + w#y + z*(x*w + y*z*(x + w)) #x*w + y*z*(x + w) + z
f1 = expand(f1)
print(f'f1: {f1}, type: {type(f1)}')
cost_original = 0
f2 = x**2+x
f3 = y**2+y
f4 = z**2+z
f5 = w**2+w

# Compute the Grobner basis
GB = groebner([f1, f2, f3, f4, f5], x, y, z, w, order='grlex')
# 'lex' (lexicographic order)
# 'grlex' (graded lexicographic order)
# 'grevlex' (graded reverse lexicographic order)
print(GB)

# add all polynomials together to make a new polynomial
new_poly = 0

#calculate the cost of f1
# for j in f1.args:
#     print(f'j: {j}, type: {type(j)}')
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
print("f1_new: ", f1)
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
# print(f1)
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
            if count<=10:
                cost_original += cost_lut[count]
                # print(f"cost_lut[count]: {cost_lut[count]}, cost_original: {cost_original}")
            else:
                cost_original += 12*count-34
elif isinstance(f1, Number):
    if f1!=0:
        cost_original += 1
        # print(f"f1 is a Number: {f1}, cost: {cost_original}")
# check the polynomials in the Grobner basis, if there is a variable with power greater than 1, then reduce the power to 1
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
    # print("new_poly is a Number")
    if new_poly!=0:
        cost += 1
elif len(new_poly.atoms())==1:
    # print(f"len(new_poly.atoms()): {len(new_poly.atoms())}, new_poly.atoms(): {new_poly.atoms()}")
    # print("len(new_poly.args[0].args): ", len(new_poly.args[0].args))
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
    print(f'modi: {modified_poly}, type: {type(modified_poly)}')
    new_poly = modified_poly
    new_poly = Poly(new_poly)
    if isinstance(new_poly, Poly):
        # print("new_poly is a Poly")
        # print(new_poly)
        if not new_poly.degree() == 0:
            # new_poly = Poly(new_poly)
            # print(new_poly)
            # for every cube in the polynomial, count the number of literals in each cube
            # print(new_poly.args[0])
            # print(len(new_poly.args[0].args))
            if len(new_poly.free_symbols) == 1 and len(new_poly.args[0].args) <= 1:
                symbol = next(iter(new_poly.free_symbols))
                # print("Symbol in the polynomial:", symbol)
                cost += 1
            for j in new_poly.args[0].args:
                count = 0
                # print(len(j.args))
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
            # print(f"new_poly is a Poly: {new_poly}, cost: {cost}")
    elif isinstance(new_poly, Number):
        if new_poly!=0:
            cost += 1
            # print(f"new_poly is a Number: {new_poly}, cost: {cost}")

print(new_poly)
print(f"cost: {cost}")
print(f"cost_original: {cost_original}")


