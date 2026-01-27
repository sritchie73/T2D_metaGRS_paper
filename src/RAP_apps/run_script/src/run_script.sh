#!/bin/bash

main() {
  
    # Set user-specified environment variables
    for uenv in ${env[@]}; do
      export $uenv
    done
  
    # Mount project storage to /mnt/project with dxfuse
    echo "Mounting project storage to /mnt/project/"
    echo "{
      \"files\" : [],
      \"directories\" : [
        {
          \"proj_id\" : \"$DX_PROJECT_CONTEXT_ID\",
          \"folder\" : \"/\",
          \"dirname\" : \"/project\"
        }
      ]
    }" > dxfuse_manifest.json
    mkdir -p /mnt
    dxfuse /mnt dxfuse_manifest.json
    rm dxfuse_manifest.json
    
    # Extract R package library bundled with app
    echo "Extracting prebundled and precompiled R packages..."
    sudo tar -xzf /usr/lib/R/site-library/Rpackages.tar.gz -C /usr/lib/R/site-library
    sudo rm -f /usr/lib/R/site-library/Rpackages.tar.gz
    
    # Install nextflow
    echo "Installing nextflow..."
    curl -s https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
    
    # Download the user supplied script and run it (or the command supplied)
    dx download "$script"
    chmod +x $script_name
    
    if [[ "$cmd" != "" ]]; then 
      echo "Running user supplied command..."
      bash -c "$cmd"
    else 
      echo "Running user supplied script..."
      ./$script_name
    fi

}
