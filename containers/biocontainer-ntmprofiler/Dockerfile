#=================================================
# PATCH EXISTING BIOCONTAINER 
#=================================================

#NOTE: Make sure the version of ntm-profiler version is correct
FROM quay.io/biocontainers/ntm-profiler:0.4.0--pyhdfd78af_0

#NOTE: Just update the ntm-profiler to the latest version of the database
RUN ntm-profiler update_db  --logging DEBUG

