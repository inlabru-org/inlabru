trd = read.csv(file="data/ETP effort.csv",sep=",")

tr = trd[1000:3000,]
tr = tr[!is.na(tr[,1]),]

plot(tr$lon, tr$lat,xlim=c(-160,0),ylim=c(-20,40))

plot(tr$lon, tr$lat)

dif = sqrt((tr$lat[1:length(tr$lat)-1] - tr$lat[2:length(tr$lat)])^2 +(tr$lon[1:length(tr$lon)-1] - tr$lon[2:length(tr$lon)])^2)
plot(dif)
hist(dif)
hist()