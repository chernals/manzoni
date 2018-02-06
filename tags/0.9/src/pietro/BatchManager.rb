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
     command_line = "bsub -J #{@job_name} -q #{@queue} #{@exe_dir}/#{"jobs/"+@job_file}"
     system(command_line)
     puts command_line
     puts "Job sent."     
   end
   
   def prepareJobFile
   file = File.new("jobs/"+@job_file, File::CREAT|File::TRUNC|File::RDWR, 0755)
   file << "export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.00/x86_64-slc5-gcc43-opt/root/lib:/afs/cern.ch/sw/lcg/contrib/gcc/4.3.5/x86_64-slc5-gcc34-opt/lib64:/afs/cern.ch/sw/lcg/contrib/mpfr/2.3.1/x86_64-slc5-gcc34-opt/lib:/afs/cern.ch/sw/lcg/contrib/gmp/4.2.2/x86_64-slc5-gcc34-opt/lib:/afs/cern.ch/user/c/chernals/public/libs/xerces-c-3.1.1/lib: \n"
   file << "cp #{@exe_dir}/#{@exe_name} . \n"
   file << "cp #{@exe_dir}/#{@sim_dir}/input.xml input.xml \n"
   file << "./#{@exe_name} -d -l2 -i input.xml\n"
   file << "tar cfz data.tar.gz *.*\n"
   file << "rm *.out\n"
   file << "rm *.pdf\n"
   file << "rm manzoni \n"
   file << "cp -R * #{@exe_dir}/#{@sim_dir} \n"
   file.close
   end
 
 
end
