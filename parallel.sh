cat strings | parallel -k -group --jobs 11 echo {} | ./jellyfish-linux count -m 15 -s 10M /dev/stdin -o /dev/stdout | ./jellyfish-linux dump /dev/stdin
