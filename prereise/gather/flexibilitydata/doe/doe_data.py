import os
import zipfile
import tarfile
import urllib.request

'''
    cleanup_directory(root)
    
Recursively cleanup a folder by deleting meaningless or empty files
'''
def cleanup_directory(root):

    all_files = os.listdir(root)
    for i in all_files:
        fp =  root+'/'+i
        if i[0] == '.' or (i[-3:]=='csv' and os.path.getsize(fp) < 0.01):
            os.remove(fp)

    all_folders = os.listdir(root)
    for i in all_folders:

        if not i[-3:] == 'csv':
            cleanup_directory(root+i)

    return

'''
    download_doe()

Download demand flexibility filters from OEDI, extract and cleanup    
'''
def download_doe():

    # create data directory 
    if not os.path.exists('data'):
        os.mkdir('data')

    # download zip data
    oedi_filter_link = "https://data.openei.org/files/180/2006weatherentireusdrfilters.tar.zip"
    urllib.request.urlretrieve(oedi_filter_link, 'data/filter.zip')

    # extract
    with zipfile.ZipFile('data/filter.zip', 'r') as fh:
        fh.extractall('data')

    # delete and further extract
    os.remove('data/filter.zip')

    for i in os.listdir('data/'):
        fh = tarfile.open('data/'+i)
        fh.extractall('data')

        fh.close()
        os.remove('data/'+i)

    # cleanup
    cleanup_directory('data/')


'''
    aggregate_doe(root, out)

Aggregate sector flexibilties by summing up the percentage flexibility from all sectors
and store to output files
'''
def aggregate_doe(root, out):
    if not os.path.exists(out):
        os.mkdir(out)

    all_folders = os.listdir(root)
    for i in all_folders:
        all_csvs = os.listdir(root+i)
        # new container for total flexibility
        file_flex = pd.read_csv(root + i + '/' + all_csvs[0], index_col=0)
        total_flex = file_flex['Flexibility'].copy()
        
        for c in all_csvs[1:]:
            fn = root + i + '/' + c
            file_flex = pd.read_csv(fn, index_col=0)
            total_flex = total_flex + file_flex['Flexibility']

        total_flex.to_csv(out_path + i + '.csv')
    
    return

