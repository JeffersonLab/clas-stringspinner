#!/usr/bin/env ruby

require 'fileutils'

# parse arguments
if ARGV.size < 4
  puts """USAGE #{$0} [NUM_EVENTS] [SEED] [OUTPUT_DIR] [clas-stringspinner_ARGS]...

  NUM EVENTS                the total number of events

  SEED                      top-level seed to randomize the job seeds;
                            0 for random seed

  OUTPUT_DIR                output directory

  clas-stringspinner_ARGS   remaining arguments are forwarded to
                            `clas-stringspinner`
  """
  exit 2
end
NumEvents       = ARGV[0].to_i
TopSeed         = ARGV[1].to_i
OutputDir       = ARGV[2]
ExeArgs         = ARGV[3..]
MaxEventsPerJob = 5000 # OSG constraint

# make sure we can find the main executable
ExeName = 'clas-stringspinner'
Exe     = File.realpath(File.join File.dirname(__FILE__), ExeName)

# calculate number of events per job
events_per_job = (NumEvents / MaxEventsPerJob).times.map{ |i| MaxEventsPerJob } + [ NumEvents % MaxEventsPerJob ]
raise "events_per_job is wrong" unless NumEvents == events_per_job.sum

# make output directory, and make sure it's empty
FileUtils.mkdir_p OutputDir
raise "directory #{OutputDir} is not empty" unless Dir.empty? OutputDir

# instantiate RNG
rng = Random.new(TopSeed==0 ? Random.new_seed : TopSeed)
PythiaMaxSeed = 900000000

# check executable arguments for forbidden arguments
["--seed", "--trig", "--num-events"].each do |forbidden_arg|
  raise "argument '#{forbidden_arg}' is forbidden" if ExeArgs.join(' ').match? forbidden_arg
end

# generate job list
File.open(File.join(OutputDir, "run.list"), 'w') do |out|
  events_per_job.each do |num_events|
    out.puts [
      Exe,
      "--num-events #{num_events}",
      "--seed #{rng.rand PythiaMaxSeed}",
      *ExeArgs
    ].join ' '
  end
end

# generate slurm script
DefaultSBatchArgs = {
  'job-name'      => 'stringspinner',
  'account'       => 'clas12',
  'partition'     => 'production',
  'time'          => '1:00:00',
  'mem-per-cpu'   => 1000,
  'ntasks'        => 1,
  'cpus-per-task' => 1,
  'output'        => "/farm_out/%u/%x_%A_%a.out",
  'error'         => "/farm_out/%u/%x_%A_%a.err",
}

