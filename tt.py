from sympy import symbols, groebner, Poly, expand

# Define the variables using an array of strings
# variables = ['x', 'y', 'z']
# Define the variables using the sympy library
x, y, z = symbols('x y z')

# define output varible name
output = 'f'


# the quantum cost
cost_lut = [0, 1, 1, 5, 13, 26, 38, 50, 62, 74, 86]
cost = 0
# make a string containing all the variables including the output variable
variables = ''
# for i in GB.args[0].args[1].args:
variables += str(x)
variables += ','
variables += str(y)
variables += ','
variables += str(z)

inputs = '.i '+variables

variables += ','
variables += str(output)

outputs = '.o '+output
constants = '.c 0'
first_line = '.v '+variables
# create a empty list of strings
list = []



# Define the polynomials
f1 = x*z+(x+1)+x*(y+1)*z #x + y*z + z + 1 #
print(f1)
f1 = expand(f1)
print(f1)
cost_original = 0
f2 = x**2+x
f3 = y**2+y
f4 = z**2+z

# Compute the Grobner basis
GB = groebner([f1, f2, f3, f4], x, y, z, order='lex')
# 'lex' (lexicographic order)
# 'grlex' (graded lexicographic order)
# 'grevlex' (graded reverse lexicographic order)
print(GB)

# add all polynomials together to make a new polynomial
new_poly = 0
# for i in GB:
#     new_poly += i
# print(new_poly)
#calculate the cost of f1
for j in f1.args:
    if j.is_Pow:
        # print(f"is_Pow: {j}")
        if j.args[1] > 1:
            # print(j.args[1])
            k = j
            for m in k.args:
                if m.is_constant():
                    k = k.replace(m, 1)
                    break
            # k.args[1] = 1
            # print('b')
            # print(k)
            # print('c')
            # print(f1.args)
            # print('d')
            f1 = f1.replace(j, k)
            # print(f1.args)
            # print('e')
            # print(f1)
            # print('f')
            # if the constant coefficient of some variable is even, then set the coeffecient to 0, else set it to 1
    # print(f"j: {j}")
    # print(type(j))
for j in f1.args:
    if j.is_constant():
        # print("j is constant")
        if j%2 == 0:
            f1 = f1.replace(j, 0)
        else:
            f1 = f1.replace(j, 1)
    elif len(j.args) >0:
        print(f"j is not constant: {j.args[0]}")
        if j.args[0] % 2 == 0:
            k = j
            print(f"k: {k}")
            # k.args[0] = 0
            for m in k.args:
                print(f"m: {m}")
                if m.is_constant():
                    k = k.replace(m, 0)
                    break
            f1 = f1.replace(j, k)
        else:
            k = j
            print(f"k: {k}")
            for m in k.args:
                print(f"m: {m}")
                if m.is_constant():
                    k = k.replace(m, 1)
                    break
            # k.args[0] = 1
            f1 = f1.replace(j, k)
if f1 != 0:
    print("f1 is not 0")
    f1 = Poly(f1)
    print(f1)
    # for every cube in the polynomial, count the number of literals in each cube
    for j in f1.args[0].args:
        line = 't'
        # print(type(j))
        print('j')
        count = 0
        temp = ''
        print(len(j.args))
        if len(j.args) == 0:
            print("len is 1")
            if(j!=1):
                temp += str(j)
                temp += ','
            
        else:
            for k in j.args:
                print(k)
                print('k')
                count += 1
                temp += str(k)
                temp += ','
                print("temp: ")
                print(temp)
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
                # print(j.args[1])
                # print('b')
                # print(j.args[0])
                k = j
                for m in k.args:
                    if m.is_constant():
                        k = k.replace(m, 1)
                        break
                print('c')
                print(i.args)
                print('d')
                i = i.replace(j, k)
                print(i.args)
                # print('e')
                # print(i)
                # print('f')
    for j in i.args:       
        # if the constant coefficient of some variable is even, then set the coeffecient to 0, else set it to 1
        if j.is_constant():
        # print("j is constant")
            if j%2 == 0:
                i = i.replace(j, 0)
            else:
                i = i.replace(j, 1)
        elif len(j.args) >0:
            print(f"j is not constant: {j.args[0]}")
            if j.args[0] % 2 == 0:
                k = j
                print(f"k: {k}")
                # k.args[0] = 0
                for m in k.args:
                    print(f"m: {m}")
                    if m.is_constant():
                        k = k.replace(m, 0)
                        break
                i = i.replace(j, k)
            else:
                k = j
                print(f"k: {k}")
                for m in k.args:
                    print(f"m: {m}")
                    if m.is_constant():
                        k = k.replace(m, 1)
                        break
            # k.args[0] = 1
            i = i.replace(j, k)
    # if i is not 0, then add a line to the list of strings that use a toffoli gate to compute the polynomial to the output variable in .tfc format
    if not i.is_constant():
        print(f"i is not 0: {i}")
        print(type(i))
        i = Poly(i)
        # for every cube in the polynomial, count the number of literals in each cube
        for j in i.args[0].args:
            line = 't'
            # print(type(j))
            # print('j')
            count = 0
            temp = ''
            # print(len(j.args))
            if len(j.args) == 0:
                # print("len is 1")
                if(j!=1):
                    temp += str(j)
                    temp += ','
                
            else:
                for k in j.args:
                    # print(k)
                    # print('k')
                    count += 1
                    temp += str(k)
                    temp += ','
                    # print("temp: ")
                    # print(temp)
            count += 1
            if count<=10:
                cost += cost_lut[count]
            else:
                cost += 12*count-34
            line += str(count)
            line += ' '
            line += temp
            line += str(output)
            list.append(line)
    elif i.is_constant() and i!=0:
        cost += 1
        line = 't1 '+str(output)
        list.append(line)
    new_poly += i

print(new_poly)
print(f"cost: {cost}")
print(f"cost_original: {cost_original}")

from sympy import symbols, Poly, Number

# Define the variables
x, y = symbols('x y')

# Define your polynomial or object
poly = Poly(5, x, y)  # Replace '5' with your polynomial

# Check if it's a Poly object and has a degree
if isinstance(poly, Poly):
    is_constant = poly.degree() == 0
else:
    print("Is it a Poly object?", False)
    if isinstance(poly, Number):
        print("It's a numeric object.")
    else:
        print("It's a different SymPy object.")
        
