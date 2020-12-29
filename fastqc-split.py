from sys import argv

filename=argv[1]
moduleFound=False

dmodule2name={'Per base sequence quality':'base-quality.txt',
        'Sequence Length':'length-dist.txt',
        'Per sequence':'read-quality.txt'}

out=''
for line in open(filename):
    if len(line)>0:
        if line[0]=='>':
            if len(out)>0:
                open(foutname,'w').write(out)
                out=''
                moduleFound=False
            for moduleName in dmodule2name.keys():
                if line.find(moduleName)>0:
                    foutname=dmodule2name[moduleName]
                    moduleFound=True
                    out=''
                    continue
        elif moduleFound:
            if line[0]!='#':
                out+=line

