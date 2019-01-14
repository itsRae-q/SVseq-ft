fin = open("HG00733_22_change_part")
fout=open("HG00733_22_change.sam","w")
for line in fin:
    fout.write(line.split()[0]+" "+line.split()[1]+" "+line.split()[2]+" "+line.split()[3]+" "+line.split()[4]+"\n"+line.split()[5]+"\n"+line.split()[6]+"\n"+line.split()[7]+"\n\n")

fin.close()
fout.close()