import os
root=r'E:\data\changyuan\免疫队列'
matches=[]
for dirpath,dirnames,filenames in os.walk(root):
    for fn in filenames:
        low=fn.lower()
        if 'supplement' in low or 'tpm' in low or 'rna' in low or 'supplementary data 2' in low or 'supplementary_data_2' in low:
            matches.append(os.path.join(dirpath,fn))
for m in matches:
    print(m)
print('Done. Found',len(matches),'matches')
