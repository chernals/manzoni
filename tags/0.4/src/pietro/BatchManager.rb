class BatchManager
   def initialize(sim_dir,input_name,queue,debug)
    puts "Creating a new batch job"
                  
    @debug = debug
    @exe_dir = Dir.pwd    
    if @debug
      @exe_name = manzonid
    else
      @exe_name = "manzoni"
    end
    @sim_dir = sim_dir
    @job_file = "job#{@sim_dir.gsub(/\//,"")}-#{Time.now.to_i}.job"
    @queue = queue
    @job_name = "Manzoni"
    @input_name = input_name
                  
    if @queue != "1nd80" && @queue != "2nd80" && @queue != "1nw80" && @queue != "2nw80" && @queue != "1nd" && @queue != "2nd" && @queue != "1nw" && @queue != "2nw"
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
   file << "cp #{@exe_dir}/#{@exe_name} . \n"
   file << "cp #{@exe_dir}/#{@sim_dir}/input.xml input.xml \n"
   file << "./#{@exe_name} -d -l2 -f 1 -i input.xml\n"
   file << "tar cfz data.tar.gz *.pdf\n"
   file << "rm *.pdf\n"
   file << "rm manzoni \n"
   file << "cp -R * #{@exe_dir}/#{@sim_dir} \n"
        file.close
   end
 
 
end
