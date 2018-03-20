#!/usr/local/rvm/rubies/ruby-2.4.1/bin/ruby
#    ## old shebang /usr/bin/env ruby
#
# This file is part of protk
# Created by Ira Cooke 14/12/2010
## SMALL EDITS BY WILLIAM BARSHOP, LABORATORY OF JAMES A WOHLSCHLEGEL (UCLA)
## 4/7/2017
#
# Runs an MS/MS search using the MSGFPlus search engine
#

$VERBOSE=nil
require 'protk/constants'
require 'protk/command_runner'
require 'protk/search_tool'
require 'protk/galaxy_stager'
require 'protk/galaxy_util'

for_galaxy = GalaxyUtil.for_galaxy?
input_stager = nil

# Setup specific command-line options for this tool. Other options are inherited from SearchTool
#
search_tool=SearchTool.new([
  :database,
  :explicit_output,
  :over_write,
  :enzyme,
  :modifications,
  :methionine_oxidation,
  :carbamidomethyl,
  :glyco,
  :acetyl_nterm,
  :instrument,
  :threads])

search_tool.jobid_prefix="p"
search_tool.option_parser.banner = "Run an MSGFPlus msms search on a set of msms spectrum input files.\n\nUsage: msgfplus_search.rb [options] file1.mzML file2.mzML ..."
search_tool.options.output_suffix="_msgfplus"

search_tool.options.enzyme=1
search_tool.options.instrument=0

# MS-GF+ doesnt support fragment tol so add this manually rather than via the SearchTool defaults
search_tool.add_value_option(:precursor_tol,"20",['-p','--precursor-ion-tol tol', 'Precursor ion mass tolerance.'])
search_tool.add_value_option(:precursor_tolu,"ppm",['--precursor-ion-tol-units tolu', 'Precursor ion mass tolerance units (ppm or Da). Default=ppm'])

search_tool.add_boolean_option(:pepxml,false,['--pepxml', 'Convert results to pepxml.'])
search_tool.add_value_option(:isotope_error_range,"0,1",['--isotope-error-range range', 'Takes into account of the error introduced by chooosing a non-monoisotopic peak for fragmentation.(Default 0,1)'])
search_tool.add_value_option(:fragment_method,0,['--fragment-method method', 'Fragment method 0: As written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD, 4: Merge spectra from the same precursor'])
search_tool.add_boolean_option(:decoy_search,false,['--decoy-search', 'Build and search a decoy database on the fly. Input db should not contain decoys if this option is used'])
search_tool.add_value_option(:protocol,0,['--protocol p', '0: NoProtocol (Default), 1: Phosphorylation, 2: iTRAQ, 3: iTRAQPhospho'])
search_tool.add_value_option(:min_pep_length,6,['--min-pep-length p', 'Minimum peptide length to consider, Default: 6'])
search_tool.add_value_option(:max_pep_length,40,['--max-pep-length p', 'Maximum peptide length to consider, Default: 40'])
search_tool.add_value_option(:min_pep_charge,2,['--min-pep-charge c', 'Minimum precursor charge to consider if charges are not specified in the spectrum file, Default: 2'])
search_tool.add_value_option(:max_pep_charge,3,['--max-pep-charge c', 'Maximum precursor charge to consider if charges are not specified in the spectrum file, Default: 3'])
search_tool.add_value_option(:num_reported_matches,1,['--num-reported-matches n', 'Number of matches per spectrum to be reported, Default: 1'])
search_tool.add_boolean_option(:add_features,false,['--add-features', 'output additional features'])
search_tool.add_value_option(:java_mem,"3500M",['--java-mem mem','Java memory limit when running the search (Default 3.5Gb)'])
search_tool.add_value_option(:enzymatic_termini,2,['--ntt n','Number of enzymatic termini (default 2)'])
search_tool.add_value_option(:num_mods,3,['--num-mods n','Number of modifications per peptide (default 3)'])
search_tool.add_value_option(:custom_aas,"",['--custom-aas ca','Custom AAs. These should be provided in a valid MSGF+ parsable string...'])

exit unless search_tool.check_options(true)

# Environment with global constants
#
genv=Constants.instance

# Set search engine specific parameters on the SearchTool object
#
#msgf_bin="#{genv.msgfplusjar}"

# We need to cope with the fact that MSGFPlus.jar might not be executable so fall back to the protk predefined path

#msgf_bin = "#{genv.msgfplus_root}/MSGFPlus.jar " if !msgf_bin
msgf_bin="/galaxy-central/tools/wohl-proteomics/MSGFcrux/MSGFPlus/MSGFPlus.jar"

throw "Could not find MSGFPlus.jar" if !msgf_bin || (msgf_bin.length==0) || !File.exist?(msgf_bin)

make_msgfdb_cmd=""

@output_suffix="_msgfplus"
@output_extension= search_tool.pepxml ? ".pep.xml" : ".mzid"

db_info=search_tool.database_info

database_path=db_info.path

database_stager=nil

if for_galaxy || Pathname.new(database_path).extname.to_s.downcase != ".fasta"
  database_stager = GalaxyUtil.stage_fasta(database_path)
  database_path = database_stager.staged_path
# # Database must have fasta extension
# if 
#   File.symlink(database_path,"#{database_path}.fasta") unless File.exists?("#{database_path}.fasta")
#   # make_msgfdb_cmd << "ln -s #{database_path} #{database_path}.fasta;"
#   database_path="#{database_path}.fasta"
end

# Database must be indexed
##But first let's go ahead and swap out and make sure we have upper case database...
system("awk \'{ if ($0 !~ />/) {print toupper($0)} else {print $0} }\' #{database_path} > #{database_path}.tmp && mv #{database_path}.tmp #{database_path}")
unless FileTest.exists?("#{database_path}.canno")
  # dbdir = Pathname.new(database_path).dirname.to_s
  tdavalue=search_tool.decoy_search ? 1 : 0;
  make_msgfdb_cmd << "java -Xmx3500M -cp #{msgf_bin} edu.ucsd.msjava.msdbsearch.BuildSA -d #{database_path} -tda #{tdavalue}; "
end

puts ARGV
puts ARGV.length
puts ARGV
throw "When --output is set only one file at a time can be run" if  ( ARGV.length> 1 ) && ( search_tool.explicit_output!=nil ) 

# Run the search engine on each input file
#
ARGV.each do |filename|

  if ( search_tool.explicit_output!=nil)
    output_path=search_tool.explicit_output
  else
    output_path=Tool.default_output_path(filename,@output_extension,search_tool.output_prefix,@output_suffix)
  end


  input_path=filename.chomp
  mzid_output_path="#{output_path}.mzid"

  if for_galaxy
    original_input_file = input_path
    original_input_path = Pathname.new("#{original_input_file}")
    random_sleep = rand()*10
    puts random_sleep.to_s
    puts "That was a sleep ^"
    sleep(random_sleep)
    input_stager = GalaxyStager.new("#{original_input_file}", :extension => '.mzML')
    input_path = input_stager.staged_path
  end



  # Only proceed if the output file is not present or we have opted to over-write it
  #
  if ( search_tool.over_write || !Pathname.new(output_path).exist? )
  
    # The basic command
    #
    cmd= "#{make_msgfdb_cmd} java -Xmx#{search_tool.java_mem} -jar #{msgf_bin} -d #{database_path} -s #{input_path} -o #{mzid_output_path} "

    #Semi tryptic peptides
    #
    cmd << " -ntt #{search_tool.enzymatic_termini}"

    #Decoy searches
    #
    tdavalue=search_tool.decoy_search ? 1 : 0;
    cmd << " -tda #{tdavalue}"

    # Precursor tolerance
    #
    cmd << " -t #{search_tool.precursor_tol}#{search_tool.precursor_tolu}"
    
    # Instrument type
    cmd << " -inst #{search_tool.instrument}"
    
    cmd << " -m #{search_tool.fragment_method}"

    cmd << " -addFeatures 1"

    cmd << " -protocol #{search_tool.protocol}"

    cmd << " -minLength #{search_tool.min_pep_length}"

    cmd << " -maxLength #{search_tool.max_pep_length}"

    cmd << " -minCharge #{search_tool.min_pep_charge}"

    cmd << " -maxCharge #{search_tool.max_pep_charge}"

    cmd << " -ti #{search_tool.isotope_error_range}"

    cmd << " -n #{search_tool.num_reported_matches}"

    # Enzyme
    #
    cmd << " -e #{search_tool.enzyme}"

    # Num Threads
    #
    cmd << " -thread #{search_tool.threads}" if search_tool.threads.to_i > 0

    mods_file_content = ""

    # Variable Modifications
    #
    p "Let's check for some variable mods!"
    p "st.fix_mods are ...:#{search_tool.fix_mods}"
    p "st.var_mods are ...:#{search_tool.var_mods}"
    p "st.custom_aas are ...:#{search_tool.custom_aas}"


    if search_tool.custom_aas !=""
        mods_file_content << "#{search_tool.custom_aas.split("___").collect { |mod| mod.lstrip.rstrip }.reject {|e| e.empty? }.join("\n")}\n"
    end

    if ( search_tool.var_mods !="" && search_tool.var_mods !~/None/) # Checking for none is to cope with galaxy input
      p "Look like we're going to run some variable mods!"
      var_mods = search_tool.var_mods.gsub(/,\d+.\d+\D*[^,]___/,"___").chomp("___").gsub("__gt__",">").gsub("___,","___").gsub(",___","___").split("___").collect { |mod| mod.lstrip.rstrip }.reject {|e| e.empty? }.join("\n")
      #var_mods << "O1,M,opt,any,Oxidation\n" if search_tool.methionine_oxidation
      #var_mods << "C2H2O,*,opt,Prot-N-term,Acetyl\n" if search_tool.acetyl_nterm
      #var_mods << "H-1N-1O1,N,opt,any,Deamidated\n" if search_tool.glyco
      var_mods = var_mods.gsub(/__lt__.*__gt__/, "")
      p "var_mods are ...:#{var_mods}"
      if ( var_mods !="" )
        mods_file_content << "#{var_mods}\n"
      end
    end

  # Fixed modifications
  #
    if ( search_tool.fix_mods !="" && search_tool.fix_mods !~/None/)
      fix_mods = search_tool.fix_mods.gsub(/,\d+.\d+\D*[^,]___/,"___").chomp("___").gsub("opt","fix").gsub("__gt__",">").gsub("___,","___").gsub(",___","___").split("___").collect { |mod| mod.lstrip.rstrip }.reject { |e| e.empty? }.join("\n")
      #fix_mods << "C2H3N1O1,C,opt,any,Carbamidomethyl\n" if search_tool.carbamidomethyl
      fix_mods = fix_mods.gsub(/__lt__.*__gt__/, "")
      if ( fix_mods !="")
        mods_file_content << "#{fix_mods}"    
      end
    end



    if ( mods_file_content != "")
      mods_path="#{output_path}.msgfplus_mods.txt"
      mods_file=File.open(mods_path,'w+')
      mods_file.write "NumMods=#{search_tool.num_mods}\n#{mods_file_content}"
      mods_file.close
      cmd << " -mod #{mods_path}"
    end
    
    # As a final part of the command we convert to pepxml
    if search_tool.pepxml
      #if search_tool.explicit_output
      cmd << ";ruby -pi.bak -e \"gsub('post=\\\"?','post=\\\"X')\" #{mzid_output_path}"
      cmd << ";ruby -pi.bak -e \"gsub('pre=\\\"?','pre=\\\"X')\" #{mzid_output_path}"
      cmd << ";ruby -pi.bak -e \"gsub('id=\\\"UnspecificCleavage\\\"','id=\\\"UnspecificCleavage\\\" name=\\\"unspecific cleavage\\\"')\" #{mzid_output_path}"
      cmd << ";idconvert #{mzid_output_path} --pepXML -o #{Pathname.new(mzid_output_path).dirname}" 

 
      pepxml_output_path = "#{mzid_output_path.chomp('.mzid')}.pepXML"

      # Fix the msms_run_summary base_name attribute
      #
      if for_galaxy
        cmd << ";ruby -pi.bak -e \"gsub(/ base_name=[^ ]+/,' base_name=\\\"#{original_input_file}\\\"')\" #{pepxml_output_path}"
      end
      #Then copy the pepxml to the final output path
      cmd << "; mv #{pepxml_output_path} #{output_path}"
    else
      cmd << "; mv #{mzid_output_path} #{output_path}"
    end

    # Up to here we've formulated the command. The rest is cleanup
    p "Running:#{cmd}"
    
    # In case the user specified background running we need to create a jobscript path
    #
    code = search_tool.run(cmd,genv)
    throw "Command failed with exit code #{code}" unless code==0

  if for_galaxy 
    input_stager.restore_references(output_path)
  end

  unless database_stager.nil?
    database_stager.restore_references(output_path)
  end

  else
    genv.log("Skipping search on existing file #{output_path}",:warn)       
  end

  # Reset this.  We only want to index the database at most once
  #
  make_msgfdb_cmd=""

end
