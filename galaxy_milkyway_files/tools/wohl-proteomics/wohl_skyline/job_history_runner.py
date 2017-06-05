#This script is responsible for crawling over Galaxy histories by taking in a 
import optparse
import json
import bioblend,sys
import multiprocessing
import pandas
import collections
from bioblend.galaxy import GalaxyInstance


parser = optparse.OptionParser()
#parser.add_option("--job_id",action="store",type="string",dest="job_id")
parser.add_option("--history_id",action="store",type="string",dest="history_id")
parser.add_option("--history_name",action="store",type="string",dest="history_name")
parser.add_option("--job_id",action="store",type="string",dest="job_id")
parser.add_option("--tool_id",action="store",type="string",dest="tool_id")
parser.add_option("--output_file",action="store",type="string",dest="output_file")
(options,args) = parser.parse_args()


galaxy_address='127.0.0.1'
galaxy_API_key='37430b18a3e4610ea243c316b293d06f'


#options.history_name=options.history_name.replace("_-_-_"," ")

gi = GalaxyInstance(galaxy_address, key=galaxy_API_key)
def convert(data):
    if isinstance(data, basestring):
        return str(data.encode('utf-8'))
    elif isinstance(data, collections.Mapping):
        return dict(map(convert, data.iteritems()))
    elif isinstance(data, collections.Iterable):
        return type(data)(map(convert, data))
    else:
        return data


def checklists(d):
    for k in list(d.keys()):
        try:
            v=d[k]
            if isinstance(v, dict) and len(v.keys())>0:
                checklists(v)
            elif len(v)==0 and not isinstance(v,str):
                #if isinstance(v, list):
                #print "trimming ",k,"...."
                del d[k]
            elif isinstance(v,dict) and len(v.keys())==0:
                #print "trimming ",k,"...."
                del d[k]
            elif isinstance(v,list) and len(v)>0:
                temp={}
                i=1
                for each in v:
                    temp["item_"+str(i)]=each
                    i+=1
                d[k]=temp
        except:
            #if k=="inputs":
            #    print k,d
            #print "Didn't find ",k,"in",d
            pass
    #return d


def findMostRecentRun(incoming_tool_id,history_id,all_jobs,interface):
    #histories=gi.histories.get_histories()
    this_history=gi.histories.show_history(history_id)
    
    #Let's make a list of all the completed datasets in this history.
    history_datasets=this_history['state_ids']['running']
    
    for this_job in all_jobs:
        if this_job['tool_id']==incoming_tool_id:
            job_details=interface.jobs.show_job(this_job['id'],full_details=True)
            for each_dataset in job_details['outputs']:
                if job_details['outputs'][each_dataset]['id'] in history_datasets:
                    return this_job['id']


def galaxyDatasetToJobMapper(dataset_id,output,jobs_client,data_client):
    if not dataset_id in output:
        current_job=jobs_client.show_job(data_client.show_dataset(dataset_id=dataset_id)['creating_job'],full_details=True)
        output[dataset_id]=current_job
    else:
        current_job=output[dataset_id]

    #try:
    #    #current_job=jobs.show_job(all_jobs_dict[dataset_id],full_details=True)
    #    current_job=all_jobs_dict[dataset_id]
    #    #print "found job",current_job['id'],"for the dataset",dataset_id
    #except:
    #    #print "We failed to find the dataset_id ",dataset_id,"during the recursive activity..."
    #    #print "it should have belonged to the job",all_jobs_dict[dataset_id]
    #    #return {}
    if (len(current_job['inputs'])==0):#base case
        return output
    for each_input_name in current_job['inputs']:
        each_input=current_job['inputs'][each_input_name]
        output.update(galaxyDatasetToJobMapper(each_input['id'],output,jobs_client,data_client))
    return output



def galaxyJobTreeEnumerator(incoming_job_id,jobs_client,data_client):#,all_jobs):
    try:
        current_job=jobs_client.show_job(incoming_job_id,full_details=True)
    except:
        print "We failed to find the job_id ",incoming_job_id,"during the preparatory steps..."
        return {}
    if (len(current_job['inputs'])==0):
        return {}
    else:
        #if all_jobs is None:
        #    all_jobs=jobs.get_jobs()
        dataset_to_job={}
        #for each_job in all_jobs:
        this_job=jobs_client.show_job(options.job_id,full_details=True)
        for each_output in this_job['outputs']:
            this_output=this_job['outputs'][each_output]
            #dataset_to_job[this_output['id']]=each_job['id']
            dataset_to_job[this_output['id']]=this_job
        #output={}
        output=dataset_to_job
        for each_input_name in current_job['inputs']:
            each_input=current_job['inputs'][each_input_name]
            output.update(galaxyDatasetToJobMapper(each_input['id'],dataset_to_job,jobs_client,data_client))
        return output
            

jobs_client=bioblend.galaxy.jobs.JobsClient(gi)
data_client=bioblend.galaxy.datasets.DatasetClient(gi)
#all_jobs=jobs.get_jobs()

#job_id=findMostRecentRun(options.tool_id,options.history_id,all_jobs,gi)
job_id=options.job_id
print "We found the Rdata node job id, which is...",job_id
dependent_jobs=galaxyJobTreeEnumerator(job_id,jobs_client,data_client)

restructured_by_tool={}

for each_job in dependent_jobs:
    this_job=dependent_jobs[each_job]
    #print this_job
    #sys.exit(2)
    if this_job['tool_id'] not in restructured_by_tool:    
        restructured_by_tool[this_job['tool_id']]={}
    #restructured_by_tool[this_job['tool_id']][this_job['id']]={} #change from "each_job" to "this_job['id']
    restructured_by_tool[this_job['tool_id']][this_job['update_time']]={} #change from "each_job" to "this_job['id']

    #restructured_by_tool[this_job['tool_id']][this_job['id']]['']
    
    if 'job_metrics' in this_job and len(this_job['job_metrics'])>0:
        restructured_by_tool[this_job['tool_id']][this_job['update_time']]['jobs_metrics']=this_job['job_metrics']
    if 'job_metrics' in this_job:
        del this_job['job_metrics'] #This information causes issues...
    if 'inputs' in this_job:
        restructured_by_tool[this_job['tool_id']][this_job['update_time']]['inputs']=this_job['inputs']
        del this_job['inputs']
    if 'outputs' in this_job:
        restructured_by_tool[this_job['tool_id']][this_job['update_time']]['outputs']=this_job['outputs']
        del this_job['outputs']
    if 'stdout' in this_job:
        restructured_by_tool[this_job['tool_id']][this_job['update_time']]['stdout']=this_job['stdout']
        del this_job['stdout']
    if 'stderr' in this_job:
        restructured_by_tool[this_job['tool_id']][this_job['update_time']]['stderr']=this_job['stderr']
        del this_job['stderr']
    if 'params' in this_job:
        restructured_by_tool[this_job['tool_id']][this_job['update_time']]['params']=this_job['params']
        del this_job['params']
    restructured_by_tool[this_job['tool_id']][this_job['update_time']]['additional info']=this_job

restructured_by_tool=convert(restructured_by_tool)
checklists(restructured_by_tool)

new_jobs_df=pandas.DataFrame.from_dict(restructured_by_tool)

new_jobs_df.to_csv(open(options.output_file,'w'),sep='\t')

with open(options.output_file,'w') as json_output_handle:
    json.dump(restructured_by_tool,json_output_handle,indent=4)
    #json.dump(dependent_jobs,json_output_handle,indent=4,separators=("\n",": "))
