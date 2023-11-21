import pynucastro as pyna
import re

rl = pyna.ReacLibLibrary()

nucleus = ["He4", "C12", "O16"]
file_name = "reaction_nuke.dat"

h_burn = rl.linking_nuclei(nucleus, with_reverse=False)
print(h_burn)

pynet = pyna.PythonNetwork(libraries=[h_burn])
print(pynet)

count_total = 0
f = open(file_name, 'w', encoding="utf-8")
for pp in pynet.rates:
    count_part = 0
    for s in pp.sets:
        count_total_out = str(count_total).zfill(6)
        f.write("[reaction_" + count_total_out + "]\n")
        f.write("rate_type  = CAS\n")
        f.write("reactants  =")
        for i in pp.reactants:
            number = re.findall("\d+\.?\d*", str(i))
            word = re.findall(r'[a-zA-Z]', str(i))
            word = ''.join(word)
            f.write(" " + str(number[0]) + str(word))
        f.write("\n")
        f.write("products   =")
        for i in pp.products:
            number = re.findall("\d+\.?\d*", str(i))
            word = re.findall(r'[a-zA-Z]', str(i))
            word = ''.join(word)
            f.write(" " + str(number[0]) + str(word))
        f.write("\n")
        energy = pp.Q * 1.602176634e-12 * 1.e6
        f.write("energy     = " + str(energy) + "\n")
        rate_out = str(s.a[0]) + ' ' + str(s.a[1]) + ' ' + str(s.a[2]) + ' ' + str(s.a[3]) + ' ' + str(s.a[4]) \
            + ' ' + str(s.a[5]) + ' ' + str(s.a[6])
        f.write("parameters = "+rate_out+"\n")
        f.write("comments   = "+str(pp)+" ( part "+str(count_part)+" )"+"\n")
        # f.write(str(pp)+'\n')
        count_total += 1
        count_part += 1
        f.write("\n")
f.close()
