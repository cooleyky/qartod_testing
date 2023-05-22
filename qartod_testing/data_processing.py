import os

def build_data_path(refdes,method,stream,prefix,folder='interim',suffix='.nc'):
    # Input: 
    #   refdes: string built from OOI site, node, and sensor for chosen dataset
    #   method: 'recovered_inst', 'recovered_host', or 'telemetered'(?) 
    #   stream: name of data stream 
    #   source: 'prod' or 'dev'
    #   folder: 'interim' (default), 'processed', 'raw', or 'external'
    #
    # Returns:
    #   ds_path: relative path to dataset from notebook folder
    
    filename = '-'.join((prefix,refdes,method,stream))+suffix              # build filename from dataset type and source

    data_folder = os.path.relpath('../data')                            # path to data folder from notebook folder

    ds_path=os.path.join(data_folder,folder,filename)                   # build full relative path 
    
    return ds_path