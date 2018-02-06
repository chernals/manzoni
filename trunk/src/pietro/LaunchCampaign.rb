#
# Pietro-Manzoni
#
# Author: Cedric Hernaslteens <cedric.hernalsteens@cern.ch>
# 
# European Organization for Nuclear Research
#
# Copyright (c) 2010+ CERN. All rights reserved.
#

require "InputProcessor.rb"
require "BatchManager.rb"
class LaunchCampaign
   def initialize(interactive_flag)
        @iflag = interactive_flag  
        # Take the arguments of the mode
        getArguments             
        # Process master input file
        processMasterInput       
        # Send the jobs
        sendJobs                  
   end
   
   def sendJobs
     job_counter = 0
     @input_processor.getPaths.each do |p|
        check_path = p + "run"
        if File.exist?(check_path)
          puts "Job already sent (#{check_path} already exists)."
        else
        puts "Sending the job for the simulation: #{p} ..."
        job = BatchManager.new(p,@file_name,@queue, @debug,@model)
        job.send
        job_counter = job_counter + 1
        end
    end
    puts "Jobs sent: " + job_counter.to_s()
   end
   
   def processMasterInput
    @input_processor = InputProcessor.new("campaigns/"+@campaign + "/" + @file_name,@campaign)
     puts "Continue ? [say 'yes']"
     is_ok = STDIN.gets.chomp
     if is_ok != "yes"
       exit
     end
    @input_processor.createDirectoriesAndInputFiles(@campaign)
    @input_processor.outputIds()
   end
   
   def getArguments 
     @debug = false
        if @iflag
            puts "Interactive mode, will prompt for each parameters..."
            puts "Campaign name:"
            @campaign = STDIN.gets.chomp
            puts "Campaign name is: " + @campaign
            puts "Name of the input file ? [input.xml]"
            @file_name = STDIN.gets.chomp
            if @file_name == ""
                @file_name = "input.xml"
            end
            puts "File name is: campaigns/" + @campaign + "/" + @file_name
            puts "Queue ?"
            @queue = STDIN.gets.chomp
            if @queue == ""
              @queue = "1nw"
            end

        else
            if(ARGV[1] == "usage" || ARGV[1] == "help")
                usage
            end
            puts "Getting the parameters from the command line arguments..."             
            @campaign = ARGV[1]
            if @campaign == nil
              usage
            end  
            puts "Campaign name is: " + @campaign
            @file_name = ARGV[2]
            if @file_name == nil
              usage
            end      
            @queue = ARGV[3]                            
            puts "File name is: campaigns/" + @campaign + "/" + @file_name
            if ARGV[4] != nil && ARGV[4] == "debug"
                @debug = true
            end
            if @debug
                puts "Debug mode is on."
            else
                puts "Debug mode is off."
            end
        end
        
   end
   
   def usage
    puts "Arguments are:"
    puts "            --- campaign : campaign name"
    puts "            --- file     : input file"
    puts "            --- queue    : queue on LxBatch"
    puts "(optional)  --- debug    : will run manzoni in debug mode"
    exit
   end
end

