'''
Created on 31/08/2015

@author: kibanez
'''

#!/usr/bin/python

import sys, re, shlex , os, string, urllib, time, math, random, subprocess, shutil

import ConfigParser

import optparse

import logging

from subprocess import Popen , PIPE



######################################################################

class OptionParser(optparse.OptionParser):

    def check_required (self, opt):

        option = self.get_option(opt)

        atrib = getattr(self.values, option.dest)
        
        if atrib is None:
#            self.error("%s option not supplied" % option)
            return False
        else:
            return True
            

######################################################################

def read_cfg_file(cfg_filename):
    
    fi = open(cfg_filename,'r')
    
    config = ConfigParser.ConfigParser()
    config.readfp(fi)
    
    hash_cfg = {}
        
    for field in config.options('INPUT'):
        hash_cfg[field] = config.get('INPUT',field)
   
    for field in config.options('OUTPUT'):
        hash_cfg[field] = config.get('OUTPUT',field)
     
    for field in config.options('SOFTWARE'):
        hash_cfg[field] = config.get('SOFTWARE',field)
        
    fi.close()
    
    return hash_cfg


#######################################################################

def run(argv=None):
    
    if argv is None: argv = sys.argv    
   
    parser = OptionParser(add_help_option=True,description="The script performs CNV estimation within the regions of interest following CoNIFER strategy")
    
    parser.add_option("--cfg",default=None,help="Config file with the complete information of the target regions and paths of the files needed for the calling",dest="f_cfg")

                    
    # Se leen las opciones aportadas por el usuario
    (options, args) = parser.parse_args(argv[1:])

    if len(argv) == 1:
        sys.exit(0)
    
    if not parser.check_required("--cfg"):
        raise IOError('The cfg file does not exist')
        
               
    try:
        
        if options.f_cfg <> None:
            
            cfg_file = options.f_cfg        
          
            if not os.path.exists(cfg_file):
                raise IOError('The file %s does not exist' % (cfg_file))
            
            hash_cfg = read_cfg_file(cfg_file)

            # INPUT           
            alignment_path  = hash_cfg.get('alignment_path','')
            l_samples  = hash_cfg.get("sample_names",'').split(',')            
            analysis_bed = hash_cfg.get('analysis_bed','')

            # OUTPUT
            rpkm_path = hash_cfg.get('rpkm_path','')
            calling_path = hash_cfg.get('calling_path','')

            # SOFTWARE
            conifer_path = hash_cfg.get('conifer_path','')
            
            if not os.path.exists(calling_path):
                os.mkdir(calling_path)
            
            if not os.path.exists(rpkm_path):
                os.mkdir(rpkm_path)
                
            if not os.path.exists(alignment_path):
                raise IOError('The alignment path does not exist. %s' % (alignment_path))

            if not os.path.isfile(conifer_path):
                raise IOError('The CoNIFER main script does not exist. %s' % (conifer_path))

            if not os.path.isfile(analysis_bed):
                raise IOError('The file does not exist. %s' % (analysis_bed))

            
                
            #Configure logger
            formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
            console = logging.StreamHandler()
            console.setFormatter(formatter)
            console.setLevel(logging.INFO)
            logger = logging.getLogger("preprocess")
            logger.setLevel(logging.INFO)
            logger.addHandler(console)
            
            l_bams = []
            for bam_f in l_samples:
                abs_path = os.path.join(alignment_path,bam_f)
                if not os.path.exists(abs_path):
                    raise IOError("The bam file does not exist. Check if the introduced path is correct: %s" %(abs_path))
                else:
                    l_bams.append(abs_path)
                
            logger.info("CNV estimation will be done in the following files: \n %s \n" %(l_bams))
            
        # 1- Compute RPKM for each file
        # python /mnt/zvol/scratch/kibanez/Conifer/conifer_v0.2.2/conifer.py  rpkm --probes /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/DE_Roche_V4_NC_refseq_s20_update_june2015.bed --input /ingemm/experiment/DisplasiaEsqueletica/2015-07-14/2015-07-22/align/1-NGS01405-MiSeq-DEv4-1-2-RocheP2-Run150715-DIS387-42870_S1_align.realign.recal.bam --output /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/RPKM/DIS387.rpkm.txt
        
        for bam_f in l_bams:
            bam_name = bam_f.split('/')[len(bam_f.split('/'))-1]
            aux = bam_name + '.rpkm.txt'
            rpkm_output = os.path.join(rpkm_path,aux)            
            #rpkm_args = 'python %s rpkm --probes %s --input %s --output %s'%(conifer_path,analysis_bed,bam_f,rpkm_output)
            rpkm_args = ['python',conifer_path,'rpkm','--probes',analysis_bed,'--input',bam_f,'--output',rpkm_output]
            logger.info(rpkm_args)
            rpkm_output = Popen(rpkm_args,stdin=PIPE,stdout=PIPE,stderr=PIPE,close_fds=True,bufsize=1)
            (trash,logdata) = rpkm_output.communicate()
            rpkm_output.wait()
             
            if logdata <> "":
                if logdata.lower().find("error") <> -1:
                    raise RuntimeError('CNV_analysis_CoNIFER.run: Error when running python rpkm \n' % (logdata))
             
        
        # 2- Computing conifer main script
        # python /mnt/zvol/scratch/kibanez/Conifer/conifer_v0.2.2/conifer.py  analyze --probes  /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/DE_Roche_V4_NC_refseq_s20_update_june2015.bed --rpkm_dir /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/RPKM --output /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/analysis.hdf5 --svd 1 --write_svals /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/singular_values.txt --plot_scree /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/screeplot.png --write_sd /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/sd_values.txt

        output_file = os.path.join(calling_path,"analysis.hdf5")
        singular_values_file = os.path.join(calling_path,"singular_values.txt")
        plot_file =  os.path.join(calling_path,"screeplot.png")
        sd_values_file =  os.path.join(calling_path,"sd_values.txt")
        svd_value = '1'
        analyze_args = ['python',conifer_path,'analyze','--probes',analysis_bed,'--rpkm_dir',rpkm_path,'--output',output_file,'--svd',svd_value,'--write_svals',singular_values_file,'--plot_scree',plot_file,'--write_sd',sd_values_file]
        #arg = 'python %s analyze --probes %s --rpkm_dir %s --output %s --svd %s --write_svals %s --plot_scree %s --write_sd %s'%(conifer_path,analysis_bed,rpkm_path,output_file,svd_value,singular_values_file,plot_file,sd_values_file)
        logger.info(analyze_args)        
        analyze_output = Popen(analyze_args,stdin=PIPE,stdout=PIPE,stderr=PIPE,close_fds=True,bufsize=1)
        (trash,logdata) = analyze_output.communicate()
        analyze_output.wait()


        if logdata <> "":
            if logdata.lower().find("error") <> -1:
                raise RuntimeError('CNV_analysis_CoNIFER.run: Error when running python analyze \n' % (logdata))


        # 3 - Calling CNVs
        # python /mnt/zvol/scratch/kibanez/Conifer/conifer_v0.2.2/conifer.py call --input /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/analysis.hdf5 --output /mnt/zvol/scratch/kibanez/samples_with_CNV/CoNIFER/DE/2015-07-15/calls.txt 

        calls_file = os.path.join(calling_path,"calls.txt")
        #arg = 'python %s call --input %s --output %s'%(conifer_path,output_file,calls_file)
        call_args = ['python',conifer_path,'call','--input',output_file,'--output',calls_file]
        logger.info(call_args)
        call_output = Popen(call_args,stdin=PIPE,stdout=PIPE,stderr=PIPE,close_fds=True,bufsize=1)
        (trash,logdata) = call_output.communicate()
        call_output.wait()

        if logdata <> "":
            if logdata.lower().find("error") <> -1:
                raise RuntimeError('CNV_analysis_CoNIFER.run: Error when running python call \n' % (logdata))

        # 4 - Creating images
        
        images_dir = os.path.join(calling_path,"images")
        if not os.path.exists(images_dir):
            os.mkdir(images_dir)
        plot_args = ['python',conifer_path,'plotcalls','--input',output_file,'--calls',calls_file,'--outputdir',images_dir]
        logger.info(plot_args)    
        plot_output = Popen(plot_args,stdin=PIPE,stdout=PIPE,stderr=PIPE,close_fds=True,bufsize=1)
        (trash,logdata) = plot_output.communicate()
        plot_output.wait()

        if logdata <> "":
            if logdata.lower().find("error") <> -1:
                raise RuntimeError('CNV_analysis_CoNIFER.run: Error when running python plotcalls \n' % (logdata))

        logger.info("CoNIFER CNV estimation done! \n")

    
        
    except:
        print >> sys.stderr , '\n%s\t%s' % (sys.exc_info()[0],sys.exc_info()[1])
        sys.exit(2)

############################################################################333

if __name__=='__main__':
    
    run()



