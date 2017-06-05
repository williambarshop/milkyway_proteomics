import sys


arguments=sys.argv[1:]

input=arguments[0]
output=arguments[1]
decoy_str=arguments[2]

writer=open(output,'wb')
with open(input,'rb')as openfile:
    for eachline in openfile:
        capture=eachline.replace("{","").replace("}","")
        split=capture.split(",")
        for each in split:
            if decoy_str in each:
                writer.write(each.strip()+"\tdecoy\n")
            else:
                writer.write(each.strip()+"\ttarget\n")
                
                
writer.close()