fin='dg1.jsp'
fout='dg1.php'
fi=open(fin,'r')
fo=open(fout,'w')
isCode=False
k=0
for ln in fi:
    k=k+1
    if '<%' in ln and '%>' in ln:
        fo.write('<?php echo "'+str(k)+'"; ?>\n')
    elif '<%' in ln:
        if isCode:
            print('problem #1 at '+str(k)+'\t'+ln)
        else:
            isCode=True
            fo.write('<?php echo "'+str(k)+'"; ?>')
    elif '%>' in ln:
        if not isCode:
            print('problem #2 at '+str(k)+'\t'+ln)
        else:
            isCode=False
    elif not isCode:
        fo.write(ln)
fi.close()
fo.close()