import os, argparse
from sympy import symbols



parser = argparse.ArgumentParser()
parser.add_argument('REAL', type=str, 
                    help='Input file (.real format)')
args = parser.parse_args()
real_file = args.REAL
circ_name = os.path.basename(real_file).split('.real')[0]
qasm_file = circ_name + '.qasm'

realfile = open(real_file,'r', encoding="utf8")

content = realfile.readlines()
content = [line.replace('\n', '') for line in content]
content = [line for line in content if not line == '']
content = [line for line in content if line[0] != '#']

constants = []
garbage = []

for info in [line for line in content if line[0] == '.'][:-2]:
    infos = [i for i in info.split(' ') if len(i) > 0]
    label = infos[0][1:]
    rest = ' '.join(infos[1:])
    print('  - {}:\t{}'.format(label, rest))

    if label == 'numvars':
      numbits = int(rest)
      
    elif label == 'variables':
      bit_dict = {}
      for idx, bit_name in enumerate(rest.split(' ')):
        bit_dict[bit_name] = idx
    elif label == 'constants':
       constants = rest
    elif label == 'garbage':
       garbage = rest

gates = [line.rstrip() for line in content if line[0] != '.']
print('  - numgates:\t{}'.format(len(gates)))

x = [symbols(('x'+str(i))) if constants[i] == '-' else int(constants[i]) for i in range(numbits)]
current_x = x
output = []
for i in range(numbits):
    if constants[i] == '-':
        output.append(x[i]**2+x[i])
print(output)

print('\nConverting...', end=' ')
for idx, g in enumerate(gates):
    g_split = g.split(' ')
    numctrl = int(g_split[0][1:]) - 1
    ctr_labels = g_split[1:-1]
    tar_label = g_split[-1]
    ctr_bits = [bit_dict[c] for c in ctr_labels]
    tar_bit = bit_dict[tar_label]
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
print(current_x)

# current output list are stored in output
# current_x stores the final values at the end of the circuit
# garbage[i] == '1' if it is garbage output, else garbage[i] == '-'