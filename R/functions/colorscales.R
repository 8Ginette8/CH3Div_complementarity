# ### =========================================================================
# ### Get the packages
# ### =========================================================================

# if(!require(wsl.plot)){
#   devtools::install_git( "http://gitlab.wsl.ch/brunp/wsl.plot_r_package.git", 
#                          credentials = git2r::cred_user_pass("username", "pw"))
#   }

# if(!require(RColorBrewer)){install.packages("RColorBrewer");library(RColorBrewer)}
# if(!require(cetcolor)){install.packages("cetcolor");library(cetcolor)}

# # Display available colors in cet package
# display_cet_all()

# # Display available colors in RColorBrewer package
# par(mfrow=c(1,1))
# display.brewer.all()

# ### =========================================================================
# ### Create some data
# ### =========================================================================

# # Define bivariate normal distribution function
# bivnorm=function(mu,sig,cor,pts){
#   prob=1/(2*pi*sig[1]*sig[2]*sqrt(1-cor^2))*
#     exp(-((pts[,1]-mu[1])^2/sig[1]^2 -
#             (2*cor*(pts[,1]-mu[1])*(pts[,2]-mu[2]))/(sig[1]*sig[2])+
#             (pts[,2]-mu[2])^2/sig[2]^2)/(2*(1-cor^2)))
#   return(prob)
# }

# # Choose mean, standard deviation and correlation of distribution
# mu=c(.5,.6)
# sig=c(.3,.1)
# cor=.6

# # Define a mesh of points for which predictions of should be made
# pts=expand.grid(seq(0,1,length.out=101),seq(0,1,length.out=101))

# # apply bivnorm function to points
# prbs=bivnorm(mu,sig,cor,pts)

# # convert points to a matrix
# mat=matrix(NA,nrow=101,ncol=101)

# for(i in 1:nrow(pts)){
#   x=round(pts[i,1]*100+1)
#   y=round(pts[i,2]*100+1)
#   mat[x,y]=prbs[i]
# }

# mat=mat/max(mat)

# ### =========================================================================
# ### Plot with linear color scales
# ### =========================================================================

# ### Make a plot
# cols=brewer.pal(8,"YlOrRd")

# par(mar=c(4,4,4,4))
# image(mat,col=cols,xaxt="n",yaxt="n")
# box()

# # Add horizontal color scale
# par(xpd=NA)
# cscl(colors=cols,
#      horiz=TRUE,
#      crds=c(0.3,0.7,1.16,1.21),
#      zrng=c(0,1),
#      at=0:4/4,
#      title="TSS",
#      lablag=.8,
#      titlag=2,
#      tickle=.2,
#      tria="l")

# # Add vertical color scale
# cscl(colors=cols,
#      horiz=FALSE,
#      crds=c(1.13,1.16,0.1,0.9),
#      zrng=c(0,1),
#      at=0:4/4,
#      labs=c("0",expression(frac(1,4)),
#             expression(frac(1,2)),
#             expression(frac(3,4)),"1"),
#      title="logit probability",
#      lablag =.3,
#      titlag=3,
#      tickle=.2,
#      cx=1,
#      tria="b")

# par(xpd=F)

### =========================================================================
### Plot with cyclic color scale
### =========================================================================

# # let's assume the surface from before is a mountain and calculate its aspect
# library(raster)
# rst=raster()
# crs(rst)="+init=epsg:2056"
# extent(rst)=extent(0,1,0,1)
# dim(rst)=dim(mat)
# values(rst)=as.vector(mat[,ncol(mat):1])
# asp=terrain(rst,"aspect")

# # Get circular color scale
# col2=cet_pal(9, name = "c2", alpha = 1)

# # Plot
# par(mar=c(4,4,4,4))
# plot(asp,col=col2,axes=FALSE,legend=FALSE,bty="n",box=FALSE)

# # Add circular color scale
# cirscl(colors=col2[1:8],
#       radius=.10,
#       center=c(.82,.18),
#       at=0:7*45,
#       zrng=c(0,360),
#       bg=TRUE,
#       tickle=.1,
#       ttl="Aspect")
