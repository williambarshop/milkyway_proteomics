import sys, os, subprocess
import pandas
from os import listdir
from os.path import isfile, join

arguments=sys.argv
arguments=arguments[1:]
saint_output=arguments[0]
new_output=arguments[1]
filter_threshold=float(arguments[2])


df=pandas.read_csv(saint_output,sep='\t')

columns=['Prey','AvgSpec','ctrlCounts']
lastcolumns=[]
for each_run in set(df['Bait']):
    columns.append(each_run+"_AvgP")
    lastcolumns.append(each_run+"_MaxP")
    lastcolumns.append(each_run+"_Replicates")

columns=columns.extend(lastcolumns)
#output_df=pandas.DataFrame(index=set(df['Prey']),columns=columns)
df.sort_values(by=['Bait'],inplace=True)
#print df.transpose()
grouped=df.groupby(['Prey'],as_index=True)


#newdf = df.groupby(['Prey']).apply(lambda x:x.set_index('Bait').to_dict('dict'))#['amount'])
#newdf.index = ['/'.join(i) for i in newdf.index]
#newdf = newdf.reset_index()
#print newdf
#sys.exit(2)
#sorted=grouped.ix(grouped['Bait'].sort().index)



lists=[]
row_labels=[]
#extras=[]
##for eachprotein in grouped.groups:#set(df['Prey']):
for name,group in grouped:
    #new_dict={}
    #prebaits= group['Bait'].tolist()
    #baits=[f.replace(" ","-") for f in prebaits]
    avgps= group['AvgP'].tolist()
    #zipped_dict=dict(zip(baits,avgps))
    zipped_dict=dict(zip(group['Bait'].tolist(),avgps))
    lists.append(zipped_dict)
    #row_labels.append(name.replace(" ","-"))
    row_labels.append(name)
    
output_df=pandas.DataFrame(lists)
output_df['Prey']=row_labels
output_df.set_index('Prey',drop=True,inplace=True)
boolean_masks=[]

filter_df=output_df.copy(deep=True)

column_itr=0
#def collapseFilter(row):
    

for eachcolumn in output_df.columns.values:
    output_df[eachcolumn]=output_df[eachcolumn].astype(float).fillna(0.0)
    boolean_masks.append((output_df[eachcolumn]>=filter_threshold))
    del filter_df[eachcolumn]
    filter_df[eachcolumn]=output_df[eachcolumn]>=filter_threshold
    #column_itr+=1
#print filter_df
print len(filter_df)

final_mask=[]
for eachprot,eachrow in output_df.iterrows():
    #print eachrow
    lastdone=True
    for eachcolumn in output_df.columns.values:
        #if not lastdone:
        #    final_mask.append(False)
        if eachrow[eachcolumn] >= filter_threshold:
            final_mask.append(True)
            lastdone=True
            break
        else:
            lastdone=False
            #print str(eachrow[eachcolumn]),"is less than",str(filter_threshold)
    if not lastdone:
        final_mask.append(False)
#print len(final_mask)
    #print any(eachrow)
output_df_filtered=output_df[final_mask]

#print boolean_masks
#print filter_df_mask

output_df_filtered.to_csv(new_output,sep='\t')
