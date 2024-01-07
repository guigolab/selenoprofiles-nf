# Use a base image with a minimal Linux distribution
FROM alpine:latest

# Set the working directory
WORKDIR /app

# Download the CLI from the provided URL
RUN wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets

# Make the downloaded file executable
RUN chmod +x datasets

