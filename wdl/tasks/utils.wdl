version 1.0

# Read sample ids and movie paths from file

import "../structs.wdl"

task read_sample_sheet {
  input {
    File sample_sheet

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = ceil(size(sample_sheet)) + 20

  command <<<
    grep -v "^#" ~{sample_sheet} | cut -f1 > sample_ids.txt
    grep -v "^#" ~{sample_sheet} | cut -f2 | sed 's/,/\t/g' > movies.tsv
    
    # are the number of rows the same for each column
    if [ "$(wc -l < sample_ids.txt)" -ne "$(wc -l < movies.tsv)" ]; then 
      echo "The number of sample id column and movies column in the sample sheet do not have the same number of rows." 1>&2
      exit 1
    fi

    # are all sample ids unique
    if [ "$(sort sample_ids.txt | uniq -d)" ]; then
      echo "The sample ids in the sample sheet are not unique." 1>&2
      exit 1
    fi

    # are all movie paths unique
    if [ "$(cat movies.tsv | sed 's/\t\t*/\n/g' | sort | uniq -d)" ]; then
      echo "The movie paths in the sample sheet are not unique." 1>&2
      exit 1
    fi
  >>>

  output {
    Array[Pair[String,Array[String]]] samples = zip(read_lines("sample_ids.txt"), read_tsv("movies.tsv"))
  }

  runtime {
    cpu: 1
    memory: "1 GB"
    disk: "~{disk_size} GB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    awsBatchRetryAttempts: runtime_attributes.max_retries
    queueArn: runtime_attributes.queue_arn
    zones: runtime_attributes.zones
    docker: "ubuntu:20.04"
  }
}
