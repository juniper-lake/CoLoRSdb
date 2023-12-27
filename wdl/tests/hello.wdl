version 1.1

workflow HelloWorld {
  call WriteGreeting
  # call WriteGreetingGPU
  # call WriteGreetingGPU2
}

task WriteGreeting {
  command {
    sleep 30; echo "Hello World"; hostname
  }
  output {
    # Write output to standard out
    File output_greeting = stdout()
  }
}

task WriteGreetingGPU {
  command {
    sleep 30; echo "Hello World"; hostname; nvidia-smi
  }
  output {
    # Write output to standard out
    File output_greeting = stdout()
  }
  runtime {
    docker: "nvidia/cuda:12.1.0-runtime-ubuntu20.04"
    gpu: true
    slurm_partition: "ml"
    slurm_constraint: "ampere"
  }
}

task WriteGreetingGPU2 {
  command {
    sleep 30; echo "Hello World"; hostname; nvidia-smi
  }
  output {
    # Write output to standard out
    File output_greeting = stdout()
  }
  runtime {
    docker: "nvidia/cuda:12.1.0-runtime-ubuntu20.04"
    gpu: true
    slurm_partition: "ml"
    gpuCount: 2
  }
}
