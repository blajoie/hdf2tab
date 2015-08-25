
from __future__ import print_function
from __future__ import division

import numpy as np
import scipy as sp
import pdb
import h5py
import sys
import argparse
import logging
import time
import gzip
import re
import os
import math
import uuid
from collections import defaultdict

verboseprint=lambda *a, **k: None

def main():

    parser=argparse.ArgumentParser(description='Extract c-data from HDF5 file into TXT (matrix.gz)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', dest='in_file', type=str, required=True, help='interaction matrix hdf5 file')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--info',dest='info', action='store_true', help='interaction matrix hdf5 file')
    parser.add_argument('-o', '--output', dest='out_file', type=str, help='interaction matrix output file')
    parser.add_argument('-wm', dest='write_mode', default='all', choices=['cis','seperate','all'], help='write mode (cis=cis only maps, seperate=all cis/trans seperate fiels, all=single matrix)')
    parser.add_argument('-z', '--zoom', dest='zoom_coords', nargs='+', type=str, default=[], help='x/y axis zoom coordinate')
    parser.add_argument('-xz', '--xzoom', dest='x_zoom_coords', nargs='+', type=str, default=[], help='x axis zoom coordinate')
    parser.add_argument('-yz', '--yzoom', dest='y_zoom_coords', nargs='+', type=str, default=[], help='y axis zoom coordinate')
    parser.add_argument('-e', '--elements', dest='bed_file', nargs='+', type=str, default=[], help='x/y axis elements bed file (bed3+ format)')
    parser.add_argument('-xe', '--xelements', dest='x_bed_file', nargs='+', type=str, default=[], help='x axis elements bed file (bed3+ format)')
    parser.add_argument('-ye', '--yelements', dest='y_bed_file', nargs='+', type=str, default=[], help='y axis elements bed file (bed3+ format)')
    parser.add_argument('--element_exten', dest='element_exten', type=int, default=0, help='bp extension for all elements, +/- bp to start/end of each element')
    parser.add_argument('-c','--chrs', dest='selected_chrs', nargs='+', type=str, default=[], help='x/y axis subset of chromosomes to extract')
    parser.add_argument('-xc','--xchrs', dest='x_selected_chrs', nargs='+', type=str, default=[], help='x axis subset of chromosomes to extract')
    parser.add_argument('-yc','--ychrs', dest='y_selected_chrs', nargs='+', type=str, default=[], help='y axis subset of chromosomes to extract')
    parser.add_argument('-b','--blocksize', dest='blocksize', type=int, default=None, help='block size of HDF5 file')
    parser.add_argument('-p', dest='precision', type=int, default=4, help='output precision (# of digits)')
    parser.add_argument('--maxdim', dest='max_dimension', type=int, default=4000 , help='maximum dimension of allxall matrix - else cis only')
    parser.add_argument('--outputchrs', dest='output_chrs',  action='store_true', help='output the chromosome list file, no matrix output')
    parser.add_argument('--outputbins', dest='output_bins',  action='store_true', help='output the bin position file, no matrix output')
    parser.add_argument('--outputfactors', dest='output_factors', action='store_true', help='output the balancing factor list file, no matrix output')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args=parser.parse_args()

    in_file=args.in_file
    verbose=args.verbose
    info=args.info
    out_file=args.out_file
    write_mode=args.write_mode
    zoom_coords=args.zoom_coords
    x_zoom_coords=args.x_zoom_coords
    y_zoom_coords=args.y_zoom_coords
    bed_file=args.bed_file
    x_bed_file=args.x_bed_file
    y_bed_file=args.y_bed_file
    element_exten=args.element_exten
    selected_chrs=args.selected_chrs
    x_selected_chrs=args.x_selected_chrs
    y_selected_chrs=args.y_selected_chrs
    blocksize=args.blocksize
    precision=args.precision
    max_dimension=args.max_dimension
    output_chrs=args.output_chrs
    output_bins=args.output_bins
    output_factors=args.output_factors
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    verbose = info if info else verbose
    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None
    format_func=("{:."+str(precision)+"f}").format
    
    verboseprint("\n",end="")
    
    in_file_name=os.path.basename(in_file)
    inhdf=h5py.File(in_file,'r')
    
    # attrs
    genome=inhdf.attrs['genome'][:]
    # datasets
    bin_positions=inhdf['bin_positions'][:]
    chr_bin_range=inhdf['chr_bin_range'][:]
    chrs=list(inhdf['chrs'][:])
    # matrix shape
    dim=inhdf['interactions'].shape
    
    # ensure symmetrical
    ensure_symmetrical(dim)
    nrow,ncol=dim
    n=nrow=ncol
    
    # calculate optimal block size
    hdf_blocksize=inhdf['interactions'].chunks[0]
    blocksize=get_blocksize(hdf_blocksize,blocksize)
    
    # default factors to 1
    factors=np.ones(n,dtype='float64')
    # get factors 
    if "balance_factors" in inhdf.keys():
        factors=inhdf['balance_factors'][:]
    else:
        verboseprint("balance factors not supplied in H5")
        verboseprint("")
        
    # get _optional_ y axis headers
    headers = np.empty(ncol)
    headers[:] = np.nan
    if "headers" in inhdf.keys():
        headers=inhdf['headers'][:]
    else:
        headers=np.array([str(i)+'|'+genome+'|'+str(chrs[bin_positions[i,0]])+':'+str(bin_positions[i,1])+'-'+str(bin_positions[i,2]) for i in np.arange(n)])
    
    # build chr lookup dict
    chr_dict={}
    for i,c in enumerate(chrs):
        chr_dict[c]=i
    
    # dump hdf information (optional)
    if info:
        dump_hdf_info(in_file,in_file_name,nrow,ncol,genome,hdf_blocksize,blocksize,chrs,chr_dict,chr_bin_range,bin_positions)
    
    if(out_file==None):
        out_file=re.sub(".hdf5", "", in_file_name)
    out_file=re.sub(".matrix", "", out_file)
    out_file=re.sub(".gz", "", out_file)

    # set up zoom chrs/dict
    x_zoom_chrs=list()
    x_zoom_dict=defaultdict(list)
    y_zoom_chrs=list()
    y_zoom_dict=defaultdict(list)
    
    # process elements file
    x_bed_file = de_dupe_list(bed_file+x_bed_file)
    y_bed_file = de_dupe_list(bed_file+y_bed_file)
    x_zoom_chrs,x_zoom_dict=subset_by_bed(x_zoom_chrs,x_zoom_dict,x_bed_file,element_exten)
    y_zoom_chrs,y_zoom_dict=subset_by_bed(y_zoom_chrs,y_zoom_dict,y_bed_file,element_exten)
    
    # process zoom coordinates
    x_zoom_coords = de_dupe_list(zoom_coords+x_zoom_coords)
    y_zoom_coords = de_dupe_list(zoom_coords+y_zoom_coords)
    x_zoom_chrs,x_zoom_dict=subset_by_coords(x_zoom_chrs,x_zoom_dict,x_zoom_coords)
    y_zoom_chrs,y_zoom_dict=subset_by_coords(y_zoom_chrs,y_zoom_dict,y_zoom_coords)
    
    # process selected chromosomes
    x_selected_chrs = de_dupe_list(selected_chrs+x_selected_chrs)
    y_selected_chrs = de_dupe_list(selected_chrs+y_selected_chrs)

    x_chrs=de_dupe_list(x_selected_chrs+x_zoom_chrs)
    y_chrs=de_dupe_list(y_selected_chrs+y_zoom_chrs)
   
    # ensure x axis chrs exist in HDF5
    x_filtered_selected_chrs=[]
    for c in x_chrs:
        verboseprint('['+c+'] not found in h5 file!') if not c in chr_dict else x_filtered_selected_chrs.append(c)
    # ensure y axis chrs exist in HDF5
    y_filtered_selected_chrs=[]
    for c in y_chrs:
        verboseprint('['+c+'] not found in h5 file!') if not c in chr_dict else y_filtered_selected_chrs.append(c)
        
    # ensure axis chrs are sorted same as H5
    x_chrs=sorted(x_filtered_selected_chrs,key=lambda i:chr_dict[i])
    y_chrs=sorted(y_filtered_selected_chrs,key=lambda i:chr_dict[i])
    
    # if no specific chrs selected, use all
    x_chrs = chrs if not len(x_chrs) else x_chrs
    y_chrs = chrs if not len(y_chrs) else y_chrs
    
    # init bin mask
    verboseprint("building bin mask ... ")
    x_bin_mask=build_bin_mask(n,x_chrs,x_zoom_dict,chr_dict,chr_bin_range,bin_positions,"x")
    y_bin_mask=build_bin_mask(n,y_chrs,y_zoom_dict,chr_dict,chr_bin_range,bin_positions,"y")
    
    verboseprint("")
    
    # if cis mode, set union x/y chr/zoom selections
    if write_mode == 'cis':
        x_chrs=y_chrs=sorted(de_dupe_list(x_chrs+y_chrs),key=lambda i:chr_dict[i])
        x_bin_mask=y_bin_mask=x_bin_mask+y_bin_mask
    
    xdim,ydim=[np.sum(x_bin_mask),np.sum(y_bin_mask)]
    if xdim*ydim > max_dimension**2:
        sys.exit('\nerror: matrix too large! %d > %d (increase --maxdim if desired)\n' % (xdim*ydim,max_dimension**2)) 
       
    if output_bins:
        write_bins(out_file,x_bin_mask,y_bin_mask,chrs,bin_positions)
    
    if output_factors:
        write_factors(out_file,x_bin_mask,y_bin_mask,chrs,bin_positions,factors)
    
    # create header list
    x_headers=headers[np.nonzero(x_bin_mask)[0]]
    y_headers=headers[np.nonzero(y_bin_mask)[0]]

    # dump hdf, all x all (one matrix)
    
    if write_mode=='all' or write_mode=='cis':
        
        verboseprint("writing tsv matrix")
        if write_mode=='all':
            verboseprint("\t",xdim,"x",ydim,sep="")
            m_out_fh=gzip.open(out_file+'.matrix.gz',"wb")
            print(str(xdim)+"x"+str(ydim)+"\t"+"\t".join(x_headers),file=m_out_fh)
      
        k=0
        for c in y_chrs:
            r=chr_bin_range[chr_dict[c]]
            
            if write_mode == 'cis':    
                tmp_y_mask=y_bin_mask[r[0]:r[1]+1]
                tmp_x_mask=np.zeros(n,dtype=bool)
                tmp_x_mask[r[0]:r[1]+1]=y_bin_mask[r[0]:r[1]+1]
            elif write_mode == 'all':
                tmp_y_mask=y_bin_mask[r[0]:r[1]+1]
                tmp_x_mask=x_bin_mask
            
            if np.sum(tmp_y_mask)==0 or np.sum(tmp_x_mask)==0:
                verboseprint("\tskipping [",c,"] - no bins overlap selection",sep="")
                continue
                
            # find min/max idx for y_axis subset mask
            subset_end=max(np.nonzero(tmp_y_mask))[-1]
            subset_start=max(np.nonzero(tmp_y_mask))[0]
            
            # increase by chr offset, translate to real genome-wide indices
            subset_start += r[0] # minimal subset - may be internal gaps
            subset_end += r[0] # maximal subset - may be internal gaps
            # always works in block increments - adjust if needed
            c_start = subset_start
            c_start -= c_start%blocksize
            c_end = subset_end
            
            if write_mode == 'cis':    
                c_ydim=np.sum(tmp_y_mask)
                c_xdim=np.sum(tmp_x_mask)
                m_out_fh=gzip.open(out_file+'__'+c+'.matrix.gz',"wb")
                verboseprint("\n\t",c_xdim,"x",c_ydim," [",c,"]",sep="")
                x_headers=headers[np.nonzero(tmp_x_mask)[0]]
                print(str(c_xdim)+"x"+str(c_ydim)+"\t"+"\t".join(x_headers),file=m_out_fh)
                
            # perform main y_axis loop
            for i in xrange(c_start,c_end+1,blocksize):
                i_offset=max(0,subset_start-i)
                # adjust ending bin if extends beyond c_end/subset_end
                b=min(c_end+1-i,blocksize-i_offset)
                
                tmp_bin_mask=np.zeros(b,dtype=bool)
                tmp_bin_mask=y_bin_mask[i+i_offset:i+i_offset+b]
                tmp_offsets=np.nonzero(tmp_bin_mask)[0]
                
                if len(tmp_offsets) == 0:
                    continue
                
                current_block=inhdf['interactions'][i+i_offset:i+i_offset+b,:][:,tmp_x_mask][tmp_bin_mask,:]
                for j in xrange(current_block.shape[0]):
                    if headers[i+i_offset+tmp_offsets[j]] != y_headers[k]:
                        sys.exit('mask error! i+tmp_offsets[j]='+str(i+i_offset+tmp_offsets[j])+' ['+headers[i+i_offset+tmp_offsets[j]]+']  k='+str(k)+' ['+y_headers[k]+']')
                    print(y_headers[k]+"\t"+"\t".join(map(format_func,current_block[j,:])),file=m_out_fh)
                    
                    pc=((float(k)/float((ydim-1)))*100)
                    verboseprint("\r",""*50,"\r\t"+str(k)+" / "+str(ydim-1)+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
                    if verbose: sys.stdout.flush()
                    k += 1
                    
            if write_mode == 'cis' or write_mode == 'seperate':
                m_out_fh.close()
        
        if write_mode == 'all':
            m_out_fh.close()
    
    elif write_mode=='seperate':
        
        verboseprint("writing tsv matrix")
        
        k=0
        for yc in y_chrs:
            yr=chr_bin_range[chr_dict[yc]]
            tmp_y_mask=y_bin_mask[yr[0]:yr[1]+1]
            
            if np.sum(tmp_y_mask)==0:
                verboseprint("\tskipping [",yc,"] - no bins overlap selection",sep="")
                continue    
            c_ydim=np.sum(tmp_y_mask)
            
            # find min/max idx for y_axis subset mask
            subset_end=max(np.nonzero(tmp_y_mask))[-1]
            subset_start=max(np.nonzero(tmp_y_mask))[0]
            
            # increase by chr offset, translate to real genome-wide indices
            subset_start += yr[0] # minimal subset - may be internal gaps
            subset_end += yr[0] # maximal subset - may be internal gaps
            # always works in block increments - adjust if needed
            c_start = subset_start
            c_start -= c_start%blocksize
            c_end = subset_end
            
            x_file_handles=dict()
            x_bin_masks=dict()
            x_offset=0
            
            verboseprint("")
            for xc in x_chrs:
                x_file_handles[xc]=gzip.open(out_file+'__'+xc+'__'+yc+'.matrix.gz',"wb")
                xr=chr_bin_range[chr_dict[xc]]
                c_xdim=np.sum(x_bin_mask[xr[0]:xr[1]+1])
                verboseprint("\t",c_xdim,"x",c_ydim," [",xc," x ",yc,"]",sep="")
                
                tmp_x_mask=np.zeros(np.sum(x_bin_mask),dtype=bool)                
                tmp_x_mask[x_offset:x_offset+c_xdim]=True
                x_bin_masks[xc]=tmp_x_mask
                tmp_x_headers=x_headers[np.nonzero(tmp_x_mask)[0]]
                print(str(c_xdim)+"x"+str(c_ydim)+"\t"+"\t".join(tmp_x_headers),file=x_file_handles[xc])
                x_offset += c_xdim
            
            # perform main y_axis loop
            for i in xrange(c_start,c_end+1,blocksize):
                # calculate i_offset (since we are forced to work in increments of blocksize)
                i_offset=max(0,subset_start-i)
                # adjust ending bin if extends beyond c_end/subset_end
                b=min(c_end+1-i,blocksize-i_offset)
                
                tmp_bin_mask=np.zeros(b,dtype=bool)
                tmp_bin_mask[0:b]=y_bin_mask[i+i_offset:i+i_offset+b]
                
                current_block=inhdf['interactions'][i+i_offset:i+i_offset+b,:][:,x_bin_mask][tmp_bin_mask,:]
                
                for j in xrange(current_block.shape[0]):
                    for xc in x_chrs:
                        if np.sum(x_bin_masks[xc])==0:
                            verboseprint("\tskipping [",xc,"] - no bins overlap selection",sep="")
                            continue 
                        if headers[i+i_offset+j] != y_headers[k]:
                            sys.exit('mask error! i+j='+str(i+i_offset+j)+' ['+headers[i+i_offset+j]+']  k='+str(k)+' ['+y_headers[k]+']')
                        print(y_headers[k]+"\t"+"\t".join(map(format_func,current_block[j,:][:,x_bin_masks[xc]])),file=x_file_handles[xc])
                        
                    pc=((float(k)/float((ydim-1)))*100)
                    verboseprint("\r",""*50,"\r\t"+str(k)+" / "+str(ydim-1)+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
                    if verbose: sys.stdout.flush()
                    k += 1
                        
            for xc in x_chrs:
                x_file_handles[xc].close()
                
    verboseprint("")
    verboseprint("")

        
def output_factors(out_file,x_bin_mask,y_bin_mask,chrs,bin_positions,factors):
    """output x/y axis ICE factors (after chr/zoom subset)
    """
    
    out_fh=open(out_file+'.xfactors',"w")
    print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
    for i in np.nonzero(x_bin_mask)[0]:
        print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2])+"\t"+str(factors[i]),file=out_fh)
    out_fh.close()
    
    out_fh=open(out_file+'.yfactors',"w")
    print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
    for i in np.nonzero(y_bin_mask)[0]:
        print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2])+"\t"+str(factors[i]),file=out_fh)
    out_fh.close()    
           
def output_bins(out_file,x_bin_mask,y_bin_mask,chrs,bin_positions):
    """output x/y axis bins (after chr/zoom subset)
    """
    
    out_fh=open(out_file+'.xbins',"w")
    print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
    for i in np.nonzero(x_bin_mask)[0]:
        print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2]),file=out_fh)
    out_fh.close()
    
    out_fh=open(out_file+'.ybins',"w")
    print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
    for i in np.nonzero(y_bin_mask)[0]:
        print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2]),file=out_fh)
    out_fh.close()    
            
def build_bin_mask(n,chrs,zoom_dict,chr_dict,chr_bin_range,bin_positions,axis=None):
    """build a 1D mask (x or y axis) based upon user chr/coor selection
    """
    
    bin_mask=np.zeros(n,dtype=bool)
   
    # build bin mask based on user chr/zoom selection
    for c in chrs:
        c_ind=chr_dict[c]
        r=chr_bin_range[chr_dict[c]]
        if c in zoom_dict:
            zoom_coord_arr=zoom_dict[c]
            for zoom_coord in zoom_coord_arr:
                tmp_bin_positions=bin_positions[r[0]:r[1]+1]
                for i,b in enumerate(tmp_bin_positions):
                    if b[2] < zoom_coord[1]: continue
                    if b[1] > zoom_coord[2]: break
                    overlap=is_overlap([zoom_coord[1],zoom_coord[2]], [b[1],b[2]])
                    if(overlap > 0):
                        bin_mask[r[0]+i]=True
        else:
            bin_mask[r[0]:r[1]+1]=True
    verboseprint("\t",axis," bin_mask\t",np.sum(bin_mask),sep="")
    
    return(bin_mask)
    
def dump_hdf_info(in_file,in_file_name,nrow,ncol,genome,hdf_blocksize,blocksize,chrs,chr_dict,chr_bin_range,bin_positions):
    """dump hdf info
    """
    
    verboseprint("inputFile",in_file,sep="\t")
    verboseprint("inputFileName",in_file_name,sep="\t")
    verboseprint("matrix shape\t",nrow," x ",ncol,sep="")
    verboseprint("assembly",genome,sep="\t")
    verboseprint("h5 chunk",hdf_blocksize,sep="\t")
    verboseprint("user chunk",blocksize,sep="\t")
    verboseprint("\nchrs",sep="\t")
    
    for i,c in enumerate(chrs):
        cbr=chr_bin_range[chr_dict[c]]
        start,end=bin_positions[cbr[0]][1],bin_positions[cbr[1]][2]
        size=(end-start)+1
        nbins=(cbr[1]-cbr[0])+1
        verboseprint("\t",i,"\t",c,":",start,"-",end,"\t(",size,")\t",cbr,"\t",nbins,sep="")
        
    verboseprint("")
    quit()
    
    
def get_blocksize(hdf_blocksize,blocksize):
    """adjust blocksize to be evenly divisible by hdf_blocksize
    """
    
    if blocksize == None:
        blocksize = hdf_blocksize
    else:
        if float(blocksize%hdf_blocksize) != 0:
            blocksize=int(math.ceil(float(blocksize)/float(hdf_blocksize)))*hdf_blocksize
            
    verboseprint("hdf_blocksize",hdf_blocksize)
    verboseprint("blocksize",blocksize)
    verboseprint("")
    
    return(blocksize)
    
def ensure_symmetrical(dim):
    """ensure nrow=ncol [symmetrical]
    """
    nrow,ncol=dim
    
    if nrow!=ncol:
        sys.exit('error: non-symmetrical matrix found!')
        
def subset_by_coords(zoom_chrs,zoom_dict,coord):
    """read UCSC coordinates, extract chr, coordinates 
    """
    
    # process zoom coordinates
    if(coord!=None):
        for z in coord:
            coord=split_coord(z)
            
            if coord==None:
                verboseprint("invalid coord",z)
                continue
                
            coord_chr,coord_start,coord_end=coord
            if coord_chr not in zoom_chrs:
                zoom_chrs += [coord_chr]
            zoom_dict[coord_chr].append(coord)
            
    return zoom_chrs,zoom_dict
        
def subset_by_bed(bed_chrs,bed_dict,bed_file,element_exten):
    """read bed file, extract chr, coordinates 
    """
    
    num_elements=0
    for b in bed_file:
        e_fh=open(b,"r") 
        for i,li in enumerate(e_fh):
            li=li.rstrip("\n")
            if li.startswith("#") or li.startswith("track"):
                continue
            
            lineList=li.split("\t")
            lineList[1]=max(1,(int(lineList[1])-element_exten))
            lineList[2]=(int(lineList[2])+element_exten)
            z=lineList[0]+':'+str(lineList[1])+'-'+str(lineList[2])
            
            bed_coord=split_coord(z)
            
            if bed_coord==None:
                verboseprint("invalid coord",z)
                continue
                
            bed_chr,bed_start,bed_end=bed_coord
            if bed_chr not in bed_chrs:
                bed_chrs += [bed_chr]
            bed_dict[bed_chr].append(bed_coord)
            num_elements += 1
        e_fh.close()
        
    return bed_chrs,bed_dict
        
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def getSmallUniqueString():  
    tmp_uniq=str(uuid.uuid4())
    tmp_uniq=tmp_uniq.split('-')[-1]
    return(tmp_uniq)
    
def bin2header(bin,genome,chrs,index=getSmallUniqueString()):
    #name|assembly|chr:start-end
    header=str(index)+'|'+genome+'|'+str(chrs[bin[0]])+':'+str(bin[1])+'-'+str(bin[2])
    return(header)

def deGroupChr(chr_id):
    return(chr_id.split('-')[0])
    
def deGroupHeader(header,extractBy="liteChr",index=getSmallUniqueString()):
    m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
    if m==None:
        sys.exit('error: incorrect input format!')
                
    bin_id,genome,chr_id,bin_start,bin_end=m.groups()
    chr_id=chr_id.split('-')[0]

    header=str(bin_id)+'|'+genome+'|'+str(chr_id)+':'+str(bin_start)+'-'+str(bin_end)
    
    return(header)
    
def split_coord(z):
    """validate and split zoom coordinate.
    coordinates must be UCSC formatted.
    e.g. chr1:500-1000
    chr(colon)start(hyphen)end where start <= end
    """
    z=z.replace(',','')
    zoom_coord=re.search(r'(\S+):(\d+)-(\d+)',z)
    
    if zoom_coord==None:
        return None
        
    zoom_chr,zoom_start,zoom_end=zoom_coord.groups()
    zoom_start=int(zoom_start)
    zoom_end=int(zoom_end)
    
    if(zoom_start > zoom_end):
        return None
        
    return [zoom_chr,zoom_start,zoom_end]
        
def de_dupe_list(input):
    """de-dupe a list, preserving order.
    """
    
    output = []
    for x in input:
        if x not in output:
            output.append(x)
    return output

def byte_to_megabyte(byte):
    """convert bytes into megabytes.
    """
    
    return round(((byte / 1000) / 1000),4) # megabyte
    # return float(((byte / 1024) / 1024),4) # mebibyte

    
def flip_intervals(a,b):
    """flip intervals, to ensure a < b
    """
    
    return(b,a)
    
def is_overlap(a, b):
    """test to for overlap between two intervals.
    """
    
    if(a[0] > a[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(a[0])+' > end '+str(a[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    if(b[0] > b[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(b[0])+' > end '+str(b[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    
    if a[0] < b[0] and a[1] > b[1]:
        return((b[1]-b[0])+1)
    
    if b[0] < a[0] and b[1] > a[1]:   
        return((a[1]-a[0])+1)
        
    if b[0] < a[0]:
        a,b=flip_intervals(a,b)
           
    return max(0, ( min(a[1],b[1]) - max(a[0],b[0]) ) ) 
    
if __name__=="__main__":
    main()

   