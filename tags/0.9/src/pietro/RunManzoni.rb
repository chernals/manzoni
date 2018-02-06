class RunManzoni   
  def initialize(interactive_flag)
    @iflag = interactive_flag        
    # Take the arguments of the mode
    getArguments
    # And run manzoni !
    runManzoni
  end
  
  def getArguments   
    @debug = false
    if @iflag
        puts "Interactive mode, will prompt for each parameters..."
        puts "Path of execution of Manzoni:"
        @path = STDIN.gets.chomp
        puts "Path is: " + @path
        puts "Name of the input file for Manzoni (relative to the path '#{@path}'):"
        @input = STDIN.gets.chomp
        puts "Input file is: " + @path + @input
        puts "Do you want to run Manzoni in debug mode ? (yes/no)"
        if STDIN.gets.chomp == "yes"
            @debug = true
        else
            @debug = false
        end
    else
        if(ARGV[1] == "usage" || ARGV[1] == "help")
            usage
        end
        puts "Getting the parameters from the command line arguments..."
        @path = ARGV[1]
        if @path == nil
            usage
        end
        puts "Path is: " + @path
        @input = ARGV[2]
        if @input == nil
            usage
        end
        puts "Input file is: " + @path + @input
        if ARGV[3] != nil && ARGV[3] == "debug"
            @debug = true
        end
        if @debug
            puts "Debug mode is on."
        else
            puts "Debug mode is off."
        end
    end
  end
  
  def runManzoni                 
    begin
      # Change the current directory
      Dir.chdir(@path) do
        # Execute manzoni
        if @debug
          command_line = "./manzoni -d -l2 -i #{@input}"
        else
          command_line = "./manzoni -l2 -i #{@input}"
        end
        system(command_line) 
      end
    rescue
      puts "Directory of Manzoni not reachable..."
      puts "Exiting."
      exit      
    end
  end
  
  def usage
    puts "Arguments are:"
    puts "            --- path  : path to the Manzoni running directory"
    puts "            --- file  : name of the input file for Manzoni"
    puts "(optional)  --- debug : will run Manzoni in debug mode"
    exit
  end
  
end
