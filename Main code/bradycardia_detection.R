library(sp)
library(plotly)
library(concaveman)
library(geometry)
library(ks)
library(R.matlab)
library(caTools)

#==== Defining all related functions ====#
source("/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/Thesis\ /R\ Code/kerneldensityestimate.R")
source("/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/Thesis\ /R\ Code/kerneldensity.R")
source("/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/Thesis\ /R\ Code/generate_event.R")
source("/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/Thesis\ /R\ Code/read_input.R")


#====== Data pre-processing : Read signal from database, generate events and peak detection ===#

#Obtain user input

input <- read_input()

#======Defining parameters ===========#

Fs <- input$freq
alpha <- input$alpha
kern <- input$kern
pan_tomp_peaks <- event_info$peakval
train_split <- input$train
test_split <- input$test
val_split <- input$validation

# Decide which MATLAB files to run

#Run MATLAB file titled "event_gen.m"
#Save the path of this file

#Generating the datasets

event_path <- "/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/My\ code/event_gen.mat"
event_info <- generate_event(event_path)



#Extract data from the .mat file
pan_tomp_peaks <- pan_tomp_peaks[2:nrow(pan_tomp_peaks),]
ann <- event_info$annotations

set.seed(42)
rows <- sample(nrow(pan_tomp_peaks))
shuffled_peaks <- pan_tomp_peaks[rows,]

#Segment to create training, testing and validation sets

train_samp <- floor(train_split * nrow(shuffled_peaks))
train_set <- shuffled_peaks[1:train_samp, ]

val_samp <- floor(val_split * nrow(shuffled_peaks))
val_set <- shuffled_peaks[(train_samp+1):(train_samp+val_samp), ]

test_samp <- floor(test_split * nrow(shuffled_peaks))
test_set <- shuffled_peaks[(train_samp+val_samp+1):nrow(shuffled_peaks),]


training_data <- as.matrix(cbind(train_set[,1], train_set[,3]))
n <- length(training_data[,1])


#Write the .csv files
write.csv(val_set, "/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/My\ code/val_set.csv")
write.csv(test_set, "/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/My\ code/test_set.csv")
#Run the cross_val_<kernel> file for the selected kernel
#Save the path to the .mat file
#Run test_labels.m
#Save the path of the .mat output

#==== leave one out cross validation ====#


h <- readMat("/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/My\ code/h_val.mat")
h_val <- h$h.final
h <- h_val[[1]]
eye <- matrix(1L,nrow = 2, ncol = 2)
eye[1,2] <- 0
eye[2,1] <- 0
H <- eye*h

neval <- 2048
tolx <- ceiling(mean(training_data[,1])/2)
xseq <- seq(min(training_data[,1])-(tolx*25),max(training_data[,1])+(tolx*25), length= neval)
toly <- ceiling(mean(training_data[,2])/2)
yseq <- seq(min(training_data[,2])-(toly*25),max(training_data[,2])+(toly*25), length= neval)

#--- Scale the data ---#
training_data[,1] <- training_data[,1]/Fs
xseq <- xseq/Fs
a <- density(training_data[,1], bw = h, kernel = kern, n = neval, from = min(xseq), to = max(xseq))
b <- density(training_data[,2], bw = h, kernel = kern, n = neval, from = min(yseq), to = max(yseq))
kerd_x <- as.matrix(a$y)
kerd_y <- as.matrix(b$y)
p_hat <- kerd_x %*% t(kerd_y)
x_points <- as.vector(a$x)
y_points <- as.vector(b$x)


#===== calculation of C_k ====#


kerd_hypothesis <- kerneldensity(training_data, h, training_data, type = kern)

est_z <- kerd_hypothesis[[4]]
h1 <- kerd_hypothesis[[1]]

z_ord <- sort(est_z, decreasing = FALSE)

# Compute k as 
k <- floor((n+1)*alpha)

# Here alpha = 0.05 since we want to predict within 95# 
confidencek <-  floor((n+1)*alpha)

# Fetch the z_(k) value from the sorted z_k set
if (k != 0){    
  z_k <- z_ord[k]} else {    
    z_k <- z_ord[1]}

org <- cbind(min(x_points),min(y_points))
K_0 <- kerneldensity(training_data, h , org, type = kern)
c_k <- z_k - ((K_0[[4]])/(n*H[1]))
c_k <- c_k[1,1]

#=====================Creating Bx===================#
#= We use three different methods to create Bx: Concave Hull, Convex Hull and alpha shape=#

#= Finding the points that > c-k =#

pts <- which(p_hat> c_k, arr.ind = T)


xpt= vector()
ypt= vector()
for (k in 1:length(pts[,1]))
{
  ypt[k]= y_points[pts[k,1]]
  xpt[k]= x_points[pts[k,2]]
}



#---- Method 2: Convex hull ---#

pol_coords <- as.matrix(data.frame(xpt,ypt))
pol_coords[,1] <- pol_coords[,1]-min(pol_coords[,1])
z <- chull(pol_coords)
coords<- pol_coords[c(z, z[1]), ]

#scaling on the obtained data ---#


coords[,1] <- Fs*coords[,1]
pol_coords[,1] <- Fs*pol_coords[,1]



#====== Testing ======#

#--- Testing the whole test set---#


brad_pts <- readMat ("/Users/sinjinim/Documents/ASU\ Semesters/Thesis\ Sem/My\ code/brad_labels.mat")
brad_pts <- brad_pts$row


test_val <- point.in.polygon(test_set[,1], test_set[,3], coords[,1], coords[,2])

#---- calculate test error ----#

lab_n <- length(test_set[,1])
lab_ind <- matrix(data = 1, nrow = lab_n, ncol = 1)
brad_pts <- as.matrix(brad_pts)

for (i in 1:length(brad_pts)) {
  val = brad_pts[i]
  lab_ind[val] = 0
  
}

test_val <- as.matrix(test_val)
error <- 0
for (j in 1:lab_n) {
  val1 <- test_val[j]
  val2 <- lab_ind[j]
  if (val1 != val2)
    error <- error+1
}

error_rate <- (error/lab_n) * 100


#-----Plotting graphs

#Plotting the KDE
ones <- matrix(1L, nrow = neval, ncol = neval)
plane_pts <- c_k * ones
p <- plot_ly(x = x_points, y = y_points, z = p_hat) %>% 
  add_surface() %>%
  add_surface(z= plane_pts, opacity=0.9)

axx <- list(title = "x values")
axy <- list(title = "y values")
axz <- list(title = "Value of estimate, p")


p <- p %>% layout(scene = list(xaxis=axx,yaxis=axy, zaxis=axz))
p

#Plotting Bx
plot(xpt*Fs,ypt, col="lightblue", xlab="x values", ylab="y values", main="Confidence set Bx for Rectangular kernel")

#Plotting convex hull of Bx
plot(pol_coords, main= "The convex hull of Bx for Rectangular kernel", col= 'lightblue', xlab="x values", ylab="y values")
lines(coords, col='black')


    
    
  




