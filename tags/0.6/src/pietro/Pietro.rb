#!/afs/cern.ch/user/c/chernals/public/langs/ruby-1.8.7/bin/ruby
require 'pathname'
$LOAD_PATH << File.dirname(Pathname.new(File.expand_path(__FILE__)).realpath.to_s).to_s
require "RunManzoni"
require "LaunchCampaign"

class Pietro
   def initialize
      puts "=============================================================================="
      puts "Pietro"
      puts "   ---    Version 0.3"
      puts "=============================================================================="
      
      # Select and execute the correct mode
      executeMode()
   end             
   
   def executeMode
      mode = ARGV[0]
      
      # Check that a mode as been provided as the first argument 
      if mode == nil
        puts "No mode selection !"
        usage()  
        exit()
      end
      
      # Switch for each modes and instanciate the corresponding class    
      interactive_flag = false
      case mode
       when "help"
         help()
       when "runManzoni"
         puts "Lets'go for the runManzoni mode..."      
         if(ARGV.size === 1) 
           interactive_flag = true
         end
         run_manzoni = RunManzoni.new(interactive_flag)
       when "launchCampaign"
         puts "Lets'go for the launchCampaign mode..."
         if(ARGV.size === 1)
           interactive_flag = true
         end
         launch_campaign = LaunchCampaign.new(interactive_flag)
       else
         puts "Invalid mode selection !"
         usage()
         exit()
      end   
   end
   
   def usage
      puts " Usage : ./pietro mode [args] [...]"
      puts "   * mode : select the mode (see bellow)"
      puts "   * args : if you don't provide args, Pietro will run in interactive mode"
      puts "=============================================================================="
      puts " Modes:"       
      puts "  * 0 ./ help"        
      puts "  * 1 ./ runManzoni"
      puts "  * 2 ./ launchCampaign"
      puts "=============================================================================="
   end
                                     
   def help
      usage()
      puts " Help"  
      puts "   * The first argument has to be one of the valid mode"
      puts "         - runManzoni is used to execute locally a single Manzoni simulation"
      puts "         - launchCampaign allows to send a simulation campaign to lxBatch"
      puts "   * If no other argument is given, Pietro will run in interactive mode"
      puts "   * In command line mode you can obtain more help"
      puts "      if the first argument is 'usage'" 
      puts "=============================================================================="
   end
end                              
                      
# Create the main Pietro object
pietro = Pietro.new
