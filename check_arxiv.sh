#!/bin/bash

ARXIV_PATH=~/buds/notes/arxiv

arxiv() {
    arxiv_files=$(ls "$ARXIV_PATH"/*.org 2>/dev/null)

    if [ -z "$arxiv_files" ]; then
        echo "$arxiv_files not found."
    else
        # Get the last file in the list
        last_file=$(echo "$arxiv_files" | awk 'END {print}')

        # Display the content of the last file using more
        emacs "$last_file" &
    fi
}

# Call the function with the provided command-line arguments
arxiv

