#
# Pietro-Manzoni
#
# Author: Cedric Hernaslteens <cedric.hernalsteens@cern.ch>
# 
# European Organization for Nuclear Research
#
# Copyright (c) 2010+ CERN. All rights reserved.
#

require "fileutils"
require "rexml/document"
require "pp"
include REXML

class InputProcessor
  def initialize(file_name,campaign)
    
    # Open the xml input file
    @campaign = campaign
    @file_name = file_name
    pp @file_name
    if !File.file?(@file_name)
        puts "Invalid file name !"
        exit
    end          
    file = File.new @file_name
    @master_doc = Document.new file 
    @id = 0
    @idList = Array.new
    # Prepare the variable searching operations
    @variables = Hash.new                                                
    findVariables()
    puts @variables.size.to_s() + " variables found"
    pp @variables
    
    # Then create the options
    createOptions()
  end
  
  def createOptions
    # Nice but cryptic ...
    # Thanks http://stackoverflow.com
    # Requires Ruby 1.9 features 'drop' and 'product' (also available in Ruby 1.8.7)
    @options = (v = @variables.map {|k, v| ([k] * v.length).zip(v) }).first.product(*v.drop(1)).map {|args| args.reduce({}) {|h, (k, v)| h.tap {|h| h[k] = v }}}
    #puts @options
  end
  
  def outputIds
    file = File.new("campaigns/"+@campaign + "/ids","a+")
    @idList.each_with_index do |item, index|
      file.write("#{index} #{item}\n")
    end
    file.close
  end
  
  def createDirectoriesAndInputFiles(campaign_name)
    @paths = Array.new
    @options.each do |option|
        path = "campaigns/#{campaign_name}/"
        file = File.new @file_name
        doc = Document.new file
        XPath.each(doc,"//*[@id = 'ID']") { |e| e.attributes['id'] = "%05d" % @id }
        @id = @id+1
        option.each do |name, value|
            # Set the path
	    if value > 1000.0
		    path += name + sprintf("%.1f",value) + "/"
	    elsif value < 0.1
		    path += name + sprintf("%1.5f",value) + "/"
	    else
            	path += name + value.to_s() + "/"
	    end
            # Modify the input file and add the id
            XPath.each(doc, "//*[. = '#{name}']") { |e| e.text = value } 
        end
        # Create the directory 
        FileUtils.mkdir_p path
        @paths << path
        @idList[@id-1] = path
        # Save the file
        file = File.new(path+"input.xml", File::CREAT|File::TRUNC|File::RDWR, 0644)
        doc.write(file)
        file.close()
    end
  end
        
  def findVariables
    @input_docs = Hash.new
    @master_doc.elements.each("manzoniSimulation/pietro/variables/var") do |el| 
      value = el.elements["value"]
      parseValues(value,el.attributes['name'])
      @input_docs[el.attributes['name']] = Array.new
    end
  end              
  
  def parseValues(value, element_name)          
    # Use regexp to parse and extract the values
    values = Array.new
    
    # Case of a list              
    if(value.text =~ /::/)    
      tmp = value.text.scan(/((-)*\d+\.\d+)/)
      tmp.each do |v|
        values << v[0].to_f()
       end  
    
    # Case of a sweep
    elsif (value.text =~ /;;/)
      start_value = value.text.scan(/^((-)*\d+\.\d+){1};;/)[0][0]
      increment_value = value.text.scan(/;;((-)*\d+\.\d+){1};;/)
      if increment_value.size == 0
        increment_value = 1                             
      else
        increment_value = increment_value[0][0]
      end
      end_value = value.text.scan(/;;((-)*\d+\.\d+){1}$/)[0][0]
      number = (end_value.to_f() - start_value.to_f()) / increment_value.to_f()    
      Integer(number+1).times do |i|
        values << start_value.to_f() + i*increment_value.to_f()
      end

    # Case of a single value
    elsif (value.text =~ /(-)*\d+\.\d+/)
      values << value.text.to_f()
      
    # Just in case             
    else
      puts "HuHu ??"
    end
        
    # Put this array in the Hash
    @variables[element_name] = values
  end
  
  def getPaths
    return @paths
  end
  
end                
