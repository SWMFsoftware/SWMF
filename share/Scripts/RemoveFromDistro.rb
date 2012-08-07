#!/usr/bin/ruby

require 'date'
require 'Open3'

NdaysWarning = 14

def colorize(text, color_code)
	  "\e[#{color_code}m#{text}\e[0m"
end

def red(text); colorize(text, 31); end
def green(text); colorize(text, 32); end

stdin, stdout, stderr, wait_thr = Open3.popen3("make -n test | fgrep Config.pl")

maketest = stdout.read
stdin.close  # stdin, stdout and stderr should be closed explicitly in this form.
stdout.close
stderr.close

#puts(" maketest = #{maketest}\n\n")

dirs = `find ./ -iname srcUser`.split("\n")
#puts("Dirs :: #{dirs}\n")
dirs.each{ |dir|
   mask = `fgrep -a "#NOTPUBLIC" #{dir}/*.f90`.split("\n")
   #puts("Mask :: #{mask}\n")
   mask.each{|line| 

	#puts("#{line}, #{dir}\n")
	ids = line.split()
	filename = ids[0].split(":").first

	# If we have expier date and email address give 
	if ids.length==1
		puts("WARNING file #{filename} has no contact persion or expiration date.\n")
        else 
        	contactname = ids[1].split(":").last
	        enddate = ids[2].split(":").last.split("/").map{|a| a.to_i}
		if Date.new(enddate[2],enddate[0],enddate[1]) <= Date.today
			filename=""
			# TODO remove NOTPUBLIC line and comit the change
	        end
		if Date.new(enddate[2],enddate[0],enddate[1])- Date.today < NdaysWarning
			daysleft = Date.new(enddate[2],enddate[0],enddate[1])- Date.today
			puts(red("File : #{filename} from #{contactname} will be public in #{daysleft} days.\n"))
		end
	end
 
	## Get with user module we are working on
	usermode = filename.split("/").last.scan(/ModUser(.*)\.f90/).flatten.first
        ## Module used in standard test, will not be removed. 
        if maketest.scan(/u=#{usermode}/).length > 0
		filename=""
		# TODO remove NOTPUBLIC line and comit the change
	end

	if filename.length >0
		File.delete(filename)
		puts(green("Removed file #{filename}.\n"))	
        end
} 
}

File.delete("./share/Scripts/RemoveFromDistro.rb")
puts(green("Removed this script :: ./share/Scripts/RemoveFromDistro.rb \n"))

