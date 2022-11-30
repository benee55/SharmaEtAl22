library(snow);library(Rmpi);library(doParallel);library(foreach);
# initialize an Rmpi environment
ns <- commandArgs(trailingOnly=TRUE)
# mpi.spawn.Rslaves(nslaves=ns)

# # send these commands to the workers
# mpi.bcast.cmd( id <- mpi.comm.rank() )
# mpi.bcast.cmd( ns <- mpi.comm.size() )
# mpi.bcast.cmd( host <- mpi.get.processor.name() )
# 
# # all workers execute this command
# mpi.remote.exec(paste("I am", id, "of", ns, "running on", host))

cl <- parallel::makeCluster(spec = ns, type="MPI")
print("Made Cluster")
doParallel::registerDoParallel(cl)
print("Registered Cluster")
outputMat<-foreach::foreach(jobNum=1:39, .combine = "c") %dopar% {
  jobNum
}
print("Parallelized Operations")
save(outputMat, file = "testNew.RData")
print("Job Ended")

# close down the Rmpi environment
# mpi.close.Rslaves(dellog = FALSE)
# mpi.exit()
