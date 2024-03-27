# from sympy import groebner, symbols
# from sympy.logic.boolalg import to_dnf, And, Or, Not, Xor
# # from sympy.abc import x0, x1, x2
# x0, x1, x2, x3, x4 = symbols('x0 x1 x2 x3 x4')
# # Define the polynomials
# # (x0)*(x1)*(x2), (x0+1)*(x1+1)*x2, x0*(x1+1)*(x2+1), (x0+1)*x1*(x2+1)
# f1 = Xor(x0, Xor(x1, x2))
# f2 = to_dnf(f1)
# # print(f2)
# # get a list of every term in the dnf
# terms = []
# for i in f2.args:
#     terms.append(i)
# # print(terms)
# # change every term to a polynomial, & into *, ~into +1
# polynomials = []
# for i in terms:
#     temp = 1
#     for j in i.args:
#         if j.is_Symbol:
#             temp *= j
#         elif j.is_Not:
#             temp *= j.args[0]+1
#     polynomials.append(temp)
# print(polynomials)
# # x1*x2 + x1*x3 + 2*x1
# # polynomials = [(x0)*(x1)*(x2), (x0+1)*(x1+1)*x2, x0*(x1+1)*(x2+1), (x0+1)*x1*(x2+1), x0**2 + x0, x1**2 + x1, x2**2 + x2]

# # combine every polynomial into one polynomial by doing so: f=f1+f2+f2*f2
# # x0 +  x1+ x2


# # add x0**2+x0, x1**2+x1, x2**2+x2 to the list of polynomials
# polynomials.append(x0**2 + x0)
# polynomials.append(x1**2 + x1)
# polynomials.append(x2**2 + x2)

# # Calculate the Groebner basis
# groebner_basis = groebner(polynomials, [x0, x1, x2], order='grlex')

# # combine every polynomial in the groebner basis into one polynomial by doing so: f=f1+f2+f2*f2
# new_poly = 0
# for i in groebner_basis:
#     product = new_poly*i
#     new_poly += i
#     new_poly += product
# new_poly = new_poly.expand()
# print(new_poly)



# # Print the Groebner basis
# print(groebner_basis)
# # x0 + x1 + x2    + (-x0  + x1 + x2)*( ( + x0 - x1 + x2)*(x0+ x1 - x2 )) + (x0 - x1 + x2)*(x0 + x1 - x2 )

# # current:  [a2 ^ a3 ^ a4, True, a2, a3, a4]
# # circuit done.
# # output_bit: a0
# # f1:  a2 ^ a3 ^ a4
# # ( a2 +  a3  + a4, a2, a3, a4, domain='ZZ')       
# # Poly(a2*a3 + a2*a4 + a2**2 + a2*a3**2 + a2*a4**2 + 2*a2 + a3**2*a4 + a3**2 + a3*a4**2 + 2*a3 + a4**2 + 2*a4, a2, a3, a4, domain='ZZ')
# # cost: 12
# # cost_original: 0
# # output_bit: a1
# # f1:  True
# # new_poly:  0
# # new_poly is a Number
# # 0
# # cost: 0
# # cost_original: 0
# # initial circuit cost:  4
# # original total cost:  0
# # minimized total cost:  12
# # a2 + a3 - a4,  + a2 - a3 + a4,   - a2 + a3 + a4

from sympy import symbols, Poly, factor, simplify

# Define your variables
x0, x1, x2, x3, x4, x5 = symbols('x0 x1 x2 x3 x4 x5')

# Create a Poly object for the given polynomial
poly_expression = Poly(x1*x2 + x1*x2*x4 + x1*x3*x5 )

# Factorize the polynomial
poly_expression = poly_expression.as_expr()

# group every two terms together
terms_group = []
for i in range(0, len(poly_expression.args), 2):
    if i == len(poly_expression.args)-1:
        terms_group.append(poly_expression.args[i])
    else:
        terms_group.append(poly_expression.args[i]+poly_expression.args[i+1])

print(terms_group)
cost_lut = [0, 1, 1, 5, 13, 29, 61, 125, 253, 509, 1021]
cost = 0
lines = []
for term in terms_group:
    factorized_term = factor(term)
    print(factorized_term)
    print(type(factorized_term))
    if factorized_term.is_Add:
        for sub_term in factorized_term.args:
            print(sub_term)
            if sub_term.is_Mul:
                for sub_sub_term in sub_term.args:
                    print(sub_sub_term)
    elif factorized_term.is_Mul:
        print(factorized_term.args[-1])
        if factorized_term.args[-1].is_Add:
            # find a variable not in the factorized_term
            new_var = ''
            for i in range(0, 6):
                if i not in factorized_term.free_symbols:
                    new_var = 'x'+str(i)
                    break
            print("new_var: ", new_var)
            
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
                    line += '\n'
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
                line += '\n'
                lines.append(line)
                if count<=19:
                        cost += cost_lut[count]
                else:
                    cost += 2**count-3

print(lines)

factored_expression = factor(poly_expression)

# Convert the factored expression back to a string
factored_string = (factored_expression)

print(factored_string)

print(type(poly_expression))
# determine if the poly_expression is of type Add
print(poly_expression.is_Add)
