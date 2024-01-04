#!/bin/bash

LOGS_PATH=~/output/logs

logs() {
    if [ -z "$1" ] || [ -z "$2" ]; then
        echo "Usage: logs <err|out> <pattern>"
        return 1
    fi

    case "$1" in
        "err")
            log_files=$(ls "$LOGS_PATH"/err.*"$2"* 2>/dev/null)
            ;;
        "out")
            log_files=$(ls "$LOGS_PATH"/out.*"$2"* 2>/dev/null)
            ;;
        "sh")
            log_files=$(ls "$LOGS_PATH"/*.sh 2>/dev/null)
            ;;
        *)
            echo "Invalid option. Use 'err', 'out' or 'sh'."
            return 1
            ;;
    esac

    if [ -z "$log_files" ]; then
        echo "$log_files not found."
    else
        # Get the last file in the list
        last_file=$(echo "$log_files" | awk 'END {print}')

        # Display the content of the last file using more
        more "$last_file"

	## Loop over all the files
	#for file in $log_files; do
        #    more "$file"
        #done
    fi
}

# Call the function with the provided command-line arguments
logs "$@"

