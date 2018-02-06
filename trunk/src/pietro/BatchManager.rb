#
# Pietro-Manzoni
#
# Author: Cedric Hernaslteens <cedric.hernalsteens@cern.ch>
# 
# European Organization for Nuclear Research
#
# Copyright (c) 2010+ CERN. All rights reserved.
#

class BatchManager
   def initialize(sim_dir,input_name,queue,debug,model)
    puts "Creating a new batch job"
                  
    @debug = debug
    @model = model
    @exe_dir = Dir.pwd    
    if @debug
      @exe_name = "manzonid"
    else
      @exe_name = "manzoni"
    end
    @sim_dir = sim_dir
    @job_file = "job#{@sim_dir.gsub(/\//,"")}-#{Time.now.to_i}.job"
    @queue = queue
    @job_name = "Manzoni"
    @input_name = input_name
                  
    if @queue != "1nh" && @queue != "8nh" && @queue != "1nd" && @queue != "2nd" && @queue != "1nw" && @queue != "2nw" && @queue != "8nm"
       puts "Invalid queue !"
       exit
     end
    prepareJobFile
   end                            
   
   def send
     puts "Sending the job to lxBatch..."
     #command_line = "bsub -J #{@job_name} -R \"(type==SLC5_64)\" -q #{@queue} #{@exe_dir}/#{"jobs/"+@job_file}"
     command_line = "bsub -G \"u_SLAP\" -J #{@job_name} -q #{@queue} #{@exe_dir}/#{"jobs/"+@job_file}"
     system(command_line)
     puts command_line
     puts "Job sent."     
   end
   
   def prepareJobFile
   file = File.new("jobs/"+@job_file, File::CREAT|File::TRUNC|File::RDWR, 0755)
   file << "cp #{@exe_dir}/#{@exe_name} . \n"
   file << "cp #{@exe_dir}/#{@sim_dir}/input.xml input.xml \n"
   file << ". /afs/cern.ch/user/c/chernals/scratch0/soft/root/bin/thisroot.sh\n"
   file << "which root-config\n"
   file << "./#{@exe_name} -d -l2 -i input.xml\n"
   file << "tar cfz data.tar.gz *\n"
   file << "rm -f *.pdf\n"
   file << "rm -f #{@exe_name} \n"
   file << "cp -R data.tar.gz #{@exe_dir}/#{@sim_dir} \n"
   file << "cp -R trapping.out #{@exe_dir}/#{@sim_dir} \n"
   file << "cp -R run #{@exe_dir}/#{@sim_dir} \n"
   file << "cp -R *.dat #{@exe_dir}/#{@sim_dir} \n"
   file.close
   end
 
 
end
