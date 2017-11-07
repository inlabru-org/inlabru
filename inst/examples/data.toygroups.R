# Load the data

data("toygroups")

# The data set is a simulation of animal groups residing in a 1D space. Their
# locations in x-space are sampled from a Cox process with intensity

ggplot(toygroups$df.intensity) + geom_line(aes(x=x,y=g.lambda))

# Adding the simulated group locations to this plot we obtain

ggplot(toygroups$df.intensity) + 
  geom_line(aes(x=x,y=g.lambda)) +
  geom_point(data = toygroups$groups, aes(x, y=0), pch="|")

# Each group has a size mark attached to it.
# These group sizes are sampled from an exponential distribution
# for which the rate parameter depends on the x-coordinate

ggplot(toygroups$groups) + 
  geom_point(aes(x= x, y = size))

ggplot(toygroups$df.rate) +
  geom_line(aes(x,rate))
