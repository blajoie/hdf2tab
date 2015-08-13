
from __future__ import print_function

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

def main():

    parser=argparse.ArgumentParser(description='Extract c-data from HDF5 file into TXT (matrix.gz)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', dest='infile', type=str, required=True, help='interaction matrix hdf5 file')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('-o', '--output', dest='outfile', type=str, help='interaction matrix output file')
    parser.add_argument('-z', '--zoom', dest='zoom_coords', nargs='+', type=str, help='zoom coordinate (can only select symmetrical subsets)')
    parser.add_argument('-e', '--elements', dest='elementsfile', type=str, help='elements bed file (bed3+ format)')
    parser.add_argument('-m','--bmem', dest='blockmem', type=int, help='block size for extracting (default=hdf chunk size)')
    parser.add_argument('-p', dest='precision', type=int, default=4, help='output precision (# of digits)')
    parser.add_argument('--info',dest='info', action='store_true', help='interaction matrix hdf5 file')
    parser.add_argument('--or', dest='output_relative', action='store_true', help='output file relative to input file path')
    parser.add_argument('--cis', dest='cis_mode', action='store_true', help='extract cis maps only')
    parser.add_argument('--chrs', dest='selected_chrs', nargs='+', default=['default'], help='subset of chromosomes to extract, [+] = all, [-] = none, zoom selected overrides --chrs')
    parser.add_argument('--elementexten', dest='elementexten', type=int, default=0, help='bp extension for all elements, +/- bp to start/end of each element')
    parser.add_argument('--maxdim', dest='max_dimension', type=int, default=4000 , help='maximum dimension of allxall matrix - else cis only')
    parser.add_argument('--outputchrs', dest='output_chrs',  action='store_true', help='output the chromosome list file, no matrix output')
    parser.add_argument('--outputbins', dest='output_bins',  action='store_true', help='output the bin position file, no matrix output')
    parser.add_argument('--outputfactors', dest='output_factors', action='store_true', help='output the balancing factor list file, no matrix output')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args=parser.parse_args()

    infile=args.infile
    info=args.info
    output_relative=args.output_relative
    verbose=args.verbose
    outfile=args.outfile
    cis_mode=args.cis_mode
    selected_chrs=args.selected_chrs
    elementexten=args.elementexten
    max_dimension=args.max_dimension
    zoom_coords=args.zoom_coords
    elementsfile=args.elementsfile
    blockmem=args.blockmem
    precision=args.precision
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
    verboseprint = print if verbose else lambda *a, **k: None
    format_func=("{:."+str(precision)+"f}").format
    
    verboseprint("\n",end="")
    
    infile_name=os.path.basename(infile)
    if output_relative:
       infile_name=infile
    inhdf=h5py.File(infile,'r')
    
    # attrs
    genome=inhdf.attrs['genome'][:]
    # datasets
    bin_positions=inhdf['bin_positions'][:]
    chr_bin_range=inhdf['chr_bin_range'][:]
    chrs=inhdf['chrs'][:]

    # matrix shape
    nrow=inhdf['interactions'].shape[0]
    ncol=inhdf['interactions'].shape[1]
    
    # calculate optimal block size
    itx_dtype=inhdf['interactions'].dtype
    itx_dtype_size=itx_dtype.itemsize
    hdf_blocksize=inhdf['interactions'].chunks[0]
    hdf_blockmem=byte_to_megabyte(hdf_blocksize*(ncol*itx_dtype_size))
    
    verboseprint("hdf_blocksize",hdf_blocksize)
    verboseprint("hdf_blockmem",hdf_blockmem)
    
    blocksize,blockmem=get_block_size(itx_dtype_size,hdf_blocksize,hdf_blockmem,blockmem,[nrow,ncol])    
    
    verboseprint("blocksize",blocksize)
    verboseprint("blockmem",blockmem)
    verboseprint("precision",precision)
    verboseprint("")
    
    # ensure symmetrical
    if nrow!=ncol:
        sys.exit('error: non-symmetrical matrix found!')
    n=nrow=ncol
    
    # default factors to 1s
    factors=np.ones(n,dtype='float64')
    # get factors 
    if "balance_factors" in inhdf.keys():
        factors=inhdf['balance_factors'][:]
    else:
        verboseprint("balance factors not supplied in H5")
        verboseprint("")
        
    # build chr lookup dict
    chr_dict={}
    for i,c in enumerate(chrs):
        chr_dict[c]=i
    
    # dump hdf info if selected
    if(info):
        verboseprint("inputFile",infile,sep="\t")
        verboseprint("inputFileName",infile_name,sep="\t")
        verboseprint("matrix shape\t",nrow," x ",ncol,sep="")
        verboseprint("assembly",genome,sep="\t")
        verboseprint("h5 chunk",inhdf['interactions'].chunks[0],sep="\t")
        verboseprint("user chunk",blocksize,sep="\t")
        verboseprint("factors",factors,sep="\t")
        
        verboseprint("\nchrs",sep="\t")
        for i,c in enumerate(chrs):
            cbr=chr_bin_range[chr_dict[c]]
            start,end=bin_positions[cbr[0]][1],bin_positions[cbr[1]][2]
            size=(end-start)+1
            nbins=(cbr[1]-cbr[0])+1
            verboseprint("\t",i,"\t",c,":",start,"-",end,"\t(",size,")\t",cbr,"\t",nbins,sep="")
            
        verboseprint("")
        quit()
    
    if(outfile==None):
        outfile=re.sub(".hdf5", "", infile_name)
    outfile=re.sub(".matrix", "", outfile)
    outfile=re.sub(".gz", "", outfile)

    subset_chrs=list()
    subset_dict=defaultdict(list)
    
    # process elements file
    num_elements=0
    if elementsfile != None:
        verboseprint("element bed file [",elementexten,"]",sep="")
        e_fh=open(elementsfile,"r") 
        
        for i,li in enumerate(e_fh):
            li=li.rstrip("\n")
            if li.startswith("#") or li.startswith("track"):
                continue
            
            lineList=li.split("\t")
            lineList[1]=max(1,(int(lineList[1])-elementexten))
            lineList[2]=(int(lineList[1])+elementexten)
            z=lineList[0]+':'+str(lineList[1])+'-'+str(lineList[2])
            
            zoom_coord=split_zoom_coord(z)
            
            if zoom_coord==None:
                verboseprint("\t",z," *invalid*",sep="")
                continue
                
            zoom_chr,zoom_start,zoom_end=zoom_coord
            subset_chrs += [zoom_chr]
            subset_dict[zoom_chr].append(zoom_coord)
            num_elements += 1
        verboseprint("\t",num_elements," elements",sep="")
        verboseprint("")
    
    # process zoom coordinates
    if(zoom_coords!=None):
        verboseprint("zoom coordinates")
        for z in zoom_coords:
            zoom_coord=split_zoom_coord(z)
            
            if zoom_coord==None:
                verboseprint("\t",z," *invalid*",sep="")
                continue
                
            zoom_chr,zoom_start,zoom_end=zoom_coord
            verboseprint("\t",zoom_chr,":",zoom_start,"-",zoom_end,sep="")
            subset_chrs += [zoom_chr]
            subset_dict[zoom_chr].append(zoom_coord)
        verboseprint("")
    
    if len(subset_dict) > 0 and len(selected_chrs) == 1 and selected_chrs[0] == 'default':
        selected_chrs=['-']
    if len(subset_dict) == 0 and len(selected_chrs) == 1 and selected_chrs[0] == 'default':
        selected_chrs=['+'] 
    
    # process selected chromosomes
    verboseprint("selected chromosomes",selected_chrs)
    if len(selected_chrs) == 1 and selected_chrs[0] == "+": # add all additional chrs
        selected_chrs=chrs
        verboseprint("\tall additional chrs")
    elif len(selected_chrs) == 1 and selected_chrs[0] == "-": # add no additional chrs
        selected_chrs=de_dupe_list(subset_chrs)
        verboseprint("\tno additional chrs")
    else: # user select mdoe
        selected_chrs=de_dupe_list(selected_chrs+subset_chrs)
    
    filtered_selected_chrs=[]
    for i,c in enumerate(selected_chrs):
        if not c in chr_dict:
           verboseprint('\tWARNING: specificed chr ['+c+'] not found in h5 file!')
        else:
            filtered_selected_chrs.append(c)

    # ensure selected chrs is sorted same as H5
    selected_chrs=sorted(filtered_selected_chrs,key=lambda x:chr_dict[x]) # ensure selected_chrs are sorted same as the HDF
   
    for c in selected_chrs:
        verboseprint("\t",c,sep="")
    
    # check to ensure at least 1 chr is selected
    if len(selected_chrs) == 0:
        sys.exit('no chromosomes selected! must select at least 1 [--chrs,--zoom]')
        
    verboseprint("")
    
    # init bin mask
    bin_mask=np.zeros(n,dtype=bool)
    # set bin mask accord to zoom/chrs
    verboseprint("building bin mask ... ")
    for c in selected_chrs:
        verboseprint("\t",c,sep="")
        c_ind=chr_dict[c]
        r=chr_bin_range[chr_dict[c]]
        
        if c in subset_dict:
            zoom_coord_arr=subset_dict[c]
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
    verboseprint("\tdone")
    verboseprint("")
    
    # warn user that (txt) matrix > max_dim row/col is _excessive_ 
    if(np.sum(bin_mask)>max_dimension):
        verboseprint("large matrix! %d > %d\n\tenforcing cis only mode!" % (np.sum(bin_mask),max_dimension))
        cis_mode=1
        
    # dump chr list file
    if output_chrs:
        verboseprint("writing chrs ...")
        out_fh=open(outfile+'.chrs',"w")
        print("# chrIndex\tchr\tucsc_coord\tchrSize\tchrBinStart\tchrBinEnd\tnBins",file=out_fh)
        for ci,c in enumerate(selected_chrs):
            c_ind=chr_dict[c]
            r=chr_bin_range[chr_dict[c]]
            if c in subset_dict:
                zoom_coord_arr=subset_dict[c]
                for zoom_coord in zoom_coord_arr:
                    tmp_bin_mask=np.zeros(n,dtype=bool)
                    tmp_bin_positions=bin_positions[r[0]:r[1]+1]
                    for bi,b in enumerate(tmp_bin_positions):
                        if b[2] < zoom_coord[1]: continue
                        if b[1] > zoom_coord[2]: break
                        overlap=is_overlap([zoom_coord[1],zoom_coord[2]], [b[1],b[2]])
                        if(overlap > 0):
                            tmp_bin_mask[r[0]+bi]=True
                    c_end=max(np.nonzero(tmp_bin_mask))[-1]
                    c_start=max(np.nonzero(tmp_bin_mask))[0]
                    start,end=bin_positions[c_start][1],bin_positions[c_end][2]
                    size=(end-start)+1
                    nbins=(c_end-c_start)+1
                    print(ci,"\t",c,"\t",c,":",start,"-",end,"\t",size,"\t",c_start,"\t",c_end,"\t",nbins,sep="",file=out_fh)
            else:
                start,end=bin_positions[r[0]][1],bin_positions[r[1]][2]
                size=(end-start)+1
                nbins=(r[1]-r[0])+1
                print(ci,"\t",c,"\t",c,":",start,"-",end,"\t",size,"\t",r[0],"\t",r[1],"\t",nbins,sep="",file=out_fh)
                
        out_fh.close()
        verboseprint("") 
        
    # dump hdf, chr x chr (multiple matrix)
    
    if(cis_mode == 1):
    
        verboseprint("cis only mode")
        
        verboseprint("")
        
        for c in selected_chrs:
            c_ind=chr_dict[c]
            r=chr_bin_range[chr_dict[c]]
            verboseprint(c,sep="")
            
            # reset bin_mask to all zeros
            bin_mask=np.zeros(n,dtype=bool)
            if c in subset_dict:
                zoom_coord_arr=subset_dict[c]
                for zoom_coord in zoom_coord_arr:
                    tmp_bin_mask=np.zeros(n,dtype=bool)
                    tmp_bin_positions=bin_positions[r[0]:r[1]+1]
                    for i,b in enumerate(tmp_bin_positions):
                        if b[2] < zoom_coord[1]: continue
                        if b[1] > zoom_coord[2]: break
                        overlap=is_overlap([zoom_coord[1],zoom_coord[2]], [b[1],b[2]])
                        if(overlap > 0):
                            bin_mask[r[0]+i]=True
                            tmp_bin_mask[r[0]+i]=True
                    c_end=max(np.nonzero(tmp_bin_mask))[-1]
                    c_start=max(np.nonzero(tmp_bin_mask))[0]
                    n2=np.sum(tmp_bin_mask[r[0]:r[1]+1])
                    #verboseprint("\tzoom subset [",c_start,"-",c_end,"] ",n2,"x",n2,sep="")
            else:
                bin_mask[r[0]:r[1]+1]=True
                n2=np.sum(bin_mask[r[0]:r[1]+1])
                verboseprint("\tchr subset [",r[0],"-",r[1]+1,"] ",n2,"x",n2,sep="")
                
            # interaction matrix output
            n2=np.sum(bin_mask)
            
            if n2 > (max_dimension*2):
                verboseprint("sub-matrix too large! [%s] %d > %d" % (c,n2,(max_dimension*2)))
                continue

            if output_bins:
                verboseprint("\twriting bins ...")
                out_fh=open(outfile+'__'+c+'.bins',"w")
                print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
                for i in np.nonzero(bin_mask)[0]:
                    print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2]),file=out_fh)
                out_fh.close()
                
            if output_factors:
                verboseprint("\twriting factors ...")
                out_fh=open(outfile+'__'+c+'.factors',"w")
                print("# binIndex\tbinChr\tbinStart\tbinEnd\ticeFactor",file=out_fh)
                for i in np.nonzero(bin_mask)[0]:
                    print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2])+"\t"+str(factors[i]),file=out_fh)
                out_fh.close()
            
            if not any([output_chrs,output_bins,output_factors]):
               
                verboseprint("\twriting (",n2,"x",n2,") matrix",sep="")
                m_out_fh=gzip.open(outfile+'__'+c+'.matrix.gz',"wb")

                header=[str(i)+'|'+genome+'|'+str(chrs[bin_positions[i,0]])+':'+str(bin_positions[i,1])+'-'+str(bin_positions[i,2]) for i in np.nonzero(bin_mask)[0]]
                print(str(n2)+"x"+str(n2)+"\t"+"\t".join(header),file=m_out_fh)
                           
                k=0
                tmp_chr_bin_mask=bin_mask[r[0]:r[1]+1]
                c_end=max(np.nonzero(tmp_chr_bin_mask))[-1]
                c_start=max(np.nonzero(tmp_chr_bin_mask))[0]
                tmp_chr_bin_mask=tmp_chr_bin_mask[c_start:c_end+1]
                c_start += r[0]
                c_end += r[0]
                            
                for i in xrange(c_start,c_end+1,blocksize):
                    b=min(c_end+1-i,blocksize)
                    tmp_bin_mask=tmp_chr_bin_mask[i-c_start:i-c_start+b]
                        
                    verboseprint("\r",""*20,"\r\tloading block (",i,":",i+b,") ... ",sep="",end="\r")
                    if verbose: sys.stdout.flush()
                    
                    current_block=inhdf['interactions'][i:i+b,:][:,bin_mask][tmp_bin_mask,:]
                    
                    for j in xrange(current_block.shape[0]):
                        print(header[k]+"\t"+"\t".join(map(format_func,current_block[j,:])),file=m_out_fh)
                        m_out_fh.flush()
                        
                        pc=((float(k)/float((n2)))*100)
                        verboseprint("\r\t"+str(k)+" / "+str(n2)+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
                        if verbose: sys.stdout.flush()
                        
                        k+=1

                m_out_fh.close()
                
                verboseprint('\r',end="")
                pc=((float(n2)/float((n2)))*100)
                verboseprint("\t"+str(n2)+" / "+str(n2)+" ["+str("{0:.2f}".format(pc))+"%] complete",end="")
                if verbose: sys.stdout.flush()
                
                verboseprint("")
                verboseprint("")
    else:
        
        # dump hdf, all x all (one matrix)
        
        verboseprint("all mode")
        
        for c in selected_chrs:
            r=chr_bin_range[chr_dict[c]]
            if c in subset_dict:
                n2=np.sum(bin_mask[r[0]:r[1]+1])
                verboseprint("\t",c," zoom subset [",r[0],"-",r[1]+1,"] ",n2,"x",n2,sep="")
            else:
                n2=np.sum(bin_mask[r[0]:r[1]+1])
                verboseprint("\t",c," chr subset [",r[0],"-",r[1]+1,"] ",n2,"x",n2,sep="")
        
        verboseprint("")
        
        # interaction matrix output
        n2=np.sum(bin_mask)
        
        if n2 > max_dimension:

            verboseprint("matrix too large! %d > %d" % (n,max_dimension))
            
        else:
           
            if output_bins:
                verboseprint("writing bins ...")
                out_fh=open(outfile+'.bins',"w")
                print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
                for i in np.nonzero(bin_mask)[0]:
                    print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2]),file=out_fh)
                out_fh.close()
                
            if output_factors:
                verboseprint("writing factors ...")
                out_fh=open(outfile+'.factors',"w")
                print("# binIndex\tbinChr\tbinStart\tbinEnd\ticeFactor",file=out_fh)
                for i in np.nonzero(bin_mask)[0]:
                    print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2])+"\t"+str(factors[i]),file=out_fh)
                out_fh.close()
            
            if not any([output_chrs,output_bins,output_factors]):
            
                verboseprint("all - writing (",n2,"x",n2,") matrix",sep="")
                m_out_fh=gzip.open(outfile+'.matrix.gz',"wb")
                
                # create header list
                header=list()
                for c in selected_chrs:
                    tmp_bin_mask=np.zeros(n,dtype=bool)
                    r=chr_bin_range[chr_dict[c]]
                    tmp_chr_bin_mask=bin_mask[r[0]:r[1]+1]
                    tmp_bin_mask[r[0]:r[1]+1]=tmp_chr_bin_mask
                    tmp_header=[str(i)+'|'+genome+'|'+str(chrs[bin_positions[i,0]])+':'+str(bin_positions[i,1])+'-'+str(bin_positions[i,2]) for i in np.nonzero(tmp_bin_mask)[0]]
                    header += tmp_header
                    
                print(str(n2)+"x"+str(n2)+"\t"+"\t".join(header),file=m_out_fh)
              
                k=0
                for c in selected_chrs:
                    r=chr_bin_range[chr_dict[c]]
                    tmp_chr_bin_mask=bin_mask[r[0]:r[1]+1]
                    c_end=max(np.nonzero(tmp_chr_bin_mask))[-1]
                    c_start=max(np.nonzero(tmp_chr_bin_mask))[0]
                    tmp_chr_bin_mask=tmp_chr_bin_mask[c_start:c_end+1]
                    c_start += r[0]
                    c_end += r[0]

                    for i in xrange(c_start,c_end+1,blocksize):
                        b=min(c_end+1-i,blocksize)
                        tmp_bin_mask=tmp_chr_bin_mask[i-c_start:i-c_start+b]
                            
                        verboseprint("\r",""*20,"\r\tloading block (",i,":",i+b,") ... ",sep="",end="\r")
                        if verbose: sys.stdout.flush()
                        
                        current_block=inhdf['interactions'][i:i+b,:][:,bin_mask][tmp_bin_mask,:]
                        
                        for j in xrange(current_block.shape[0]):
                            print(header[k]+"\t"+"\t".join(map(format_func,current_block[j,:])),file=m_out_fh)
                            m_out_fh.flush()
                            
                            pc=((float(k)/float((n2)))*100)
                            verboseprint("\r\t"+str(k)+" / "+str(n2)+" ["+str("{0:.2f}".format(pc))+"%] complete ... ",end="\r")
                            if verbose: sys.stdout.flush()
                            
                            k+=1
                    
                m_out_fh.close()
                
                verboseprint('\r',end="")
                pc=((float(n2)/float((n2)))*100)
                verboseprint("\t"+str(n2)+" / "+str(n2)+" ["+str("{0:.2f}".format(pc))+"%] complete",end="")
                if verbose: sys.stdout.flush()
                
                verboseprint("")
                verboseprint("")

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
    
def split_zoom_coord(z):
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
        
def get_block_size(itx_dtype_size,hdf_blocksize,hdf_blockmem,blockmem,hdf_shape):
    """choose optimal hdf block size, based on user blockmem request.
    pixel dtype * ncol = mem of single row.
    Adjust block size (num of rows) to acheive blockmem mem usage
    """
    nrow,ncol=hdf_shape
    
    if blockmem==None:
        blocksize=hdf_blocksize
        blockmem=byte_to_megabyte(blocksize*(ncol*itx_dtype_size))
    else:
        blockmem=blockmem/2
        blockmem -= 256 # reserve 256MB for non-matrix items
        
        num_block_chunks=math.floor(blockmem/hdf_blockmem)
        if(num_block_chunks > 1):
            blocksize = int(hdf_blocksize*num_block_chunks)
        else:
            blocksize = hdf_blocksize
        blockmem=byte_to_megabyte(blocksize*(ncol*itx_dtype_size))
    return blocksize,blockmem

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
    # return round(((byte / 1024) / 1024),4) # mebibyte

    
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

   