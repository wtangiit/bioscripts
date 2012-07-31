#!/usr/bin/env python
'''run fgs test while recording necessary info for the specific test case''' 

import os
import sys
from optparse import OptionParser
import ConfigParser
import subprocess
import datetime

rec_fields = ['date', 'path_fgs', 'git_tag', 'path_train', 'path_input', 'path_output']

if __name__ == '__main__':
    usage  = '''usage: %prog -c <config file>'''
    parser = OptionParser(usage)
    parser.add_option("-c", "--config",   dest="config_file", default=None, help="config file")
    parser.add_option("-t", "--tag", dest="tag", default=None, help="recompile the specified tag")
    parser.add_option("-q", "--quicktest", dest="quicktest", action="store_true", default=False,  help="quick test")
    
    (opts, args) = parser.parse_args()
    if not (opts.config_file and os.path.isfile(opts.config_file) ):
        parser.error("Missing config file: %s"%(opts.config_file, ))
        sys.exit(1)
        
    config = ConfigParser.ConfigParser()
    config.read(opts.config_file)
    
    path_fgs = config.get('program', 'path_fgs')
    path_train = config.get('training', 'path_train')
    path_input = config.get('testing', 'path_input')
    path_output = config.get('testing', 'path_output')
    git_tag = config.get('program', 'git_tag')
    
    #if not git_tag specified, checkout master
    if git_tag=="":
        git_tag = "master"
    
    print "path_fgs:", path_fgs
    print "git_tag:", git_tag
    print "path_train:", path_train
    print "path_input:", path_input
    print "path_output:", path_output
    
    rec_dict = {}
    
    rec_dict['path_fgs'] = path_fgs
    rec_dict['git_tag'] = git_tag
    rec_dict['path_train'] = path_train
    rec_dict['path_input'] = path_input
    
    rec_dict['path_output'] = path_output
    
        
#rebuild fgs exectable
    os.chdir(path_fgs)
    subprocess.call(['pwd'])
    subprocess.call(['ls', '-l'])
    
    if opts.tag:
        subprocess.call(['git', 'checkout', '%s' % opts.tag])
    else:
        subprocess.call(['git', 'checkout', '%s' % git_tag])
    
    subprocess.call(['make', 'clean'])
    subprocess.call(['make', 'fgs'])
    
#re-direct training data
    subprocess.call(['rm', 'train'])
    subprocess.call(['ln', '-s', '%s' % path_train, 'train'])
    #ln -s path_train train
    
#run test
    try:
        os.stat(path_output)
    except:
        os.makedirs(path_output)
        
    if opts.quicktest:
        subprocess.call(['tesFGS.pl', '--input=%s' % path_input,  '--stem=%s' % path_output, '--quicktest'])
        print "testFGS.pl --input=%s --stem=%s --quicktest"  %  (path_input, path_output) 
    else:
        subprocess.call(['tesFGS.pl', '--input=%s' % path_input,  '--stem=%s' % path_output])
        print "testFGS.pl --input=%s --stem=%s"  %  (path_input, path_output)
      
#record metadata

    os.chdir(path_output)
    
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    rec_dict['date'] = now
    
    rec_file = open("test.rec", 'w')
    
    for field in rec_fields:
        value = rec_dict[field]
        line = "%s: %s\n" % (field, value)
        rec_file.write(line)
        
    rec_file.close()        
    
