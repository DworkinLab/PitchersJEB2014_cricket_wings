
## ----load_libraries_&_read_in_functions, results='hide', echo=FALSE, message=FALSE----
##############################
# Functions and analyses for Pitchers et al. 2014 "The potential influence of morphology on the evolutionary divergence of an acoustic signal"

rm( list= ls())

require( car )
require( MASS )
require( boot )
require( geomorph )
require( abind )
require( MCMCglmm )
require( gdata )
require( plyr )
require( effects )
require( combinat )

#############################
##### read in functions #####
#############################

# 'apply'-able HPDinterval function
hpd <- function( X ) {
  HPDinterval( as.mcmc( X ), probs=c(0.025, 0.975) )[1:2]
}

# univariate Rsquared and partial Rsquared functions
Rsq <- function( model ){
	fitted.variance <- var(model$fitted)
	total.variance	<- var(model$fitted) + var(model$resid)
	fitted.variance / total.variance
}


PRsq <- function( model ){
	residual.variance <- var(model$resid)
	variables <- attr(terms(model), "term.labels")
		model.length <- length(variables)
		variable.name <- rep(NA, model.length )
		partial.Rsq <- rep(NA, model.length )
		univariate.Model.Rsq <- rep(NA, model.length )
			
	for (i in 1:model.length){
		variable.name[i] <- variables[i]
		drop <- parse( text=variables[i] )
		new.formula <- as.formula( paste( ".~.-", variables[i], sep=""))
		new.model <- update(model, new.formula )
		partial.Rsq[i] <- (var(new.model$resid) - residual.variance)/ var(new.model$resid)
		
		new.formula.univariate <- as.formula( paste( ".~", variables[i], sep=""))
		univariate.model <- update(model, new.formula.univariate)
		univariate.Model.Rsq[i] <- summary(univariate.model)$r.sq
		}
	
	R2 <- Rsq( model )
	adj.R2 <- summary(model)$adj.r
	
	partials <- data.frame(partial.Rsq, univariate.Model.Rsq )
	row.names(partials) <- variable.name
	
	list(FullModelRsquared=R2, FullModelAdjustedR2 = adj.R2, partials=partials	)
}



# multivariate Rsquared and partial Rsquared calculators

shapeRsq <- function( model ){
	fitted.variance <- sum(diag(var(model$fitted)))
	total.variance	<- sum(diag(var(model$fitted + model$resid)))
	fitted.variance / total.variance
}


shapePRsq <- function( model ){
# Based on the derivation from page 269 of Kutner et. al. (Applied linear statistical models edition 5.)
	residual.variance <- var(model$resid)
	variables <- attr(terms(model), "term.labels")
		model.length <- length(variables)
		variable.name <- rep(NA, model.length )
		partial.Rsq <- rep(NA, model.length )
			
	for (i in 1:model.length){
		variable.name[i] <- variables[i]
		drop <- parse( text=variables[i] )
		new.formula <- as.formula( paste( ".~.-", variables[i], sep=""))
		new.model <- update(model, new.formula )
		partial.Rsq[i] <- (sum ( diag( var(new.model$resid))) - sum( diag( residual.variance)) ) / sum( diag( var(new.model$resid)))
		}
	R2 <- shapeRsq( model )
	list(Rsquared=R2, partials=data.frame( cbind( variable.name, partial.Rsq ))	)
}

# calculate tangent approximaton for tangent approximates Procrustes Distance (Euclidean Distance)
PD <- function(x) { 
	sqrt(t(x)%*%x)}

####### When the vectors can be of arbitrary sign, use this which computes the magnitude of the vector correlation, and then computes the angle.

ang.vec.abs <- function(vec1, vec2){
	vec.cor <- abs((t(vec1) %*% vec2)/(PD(vec1)*PD(vec2)))
	vec.angle <- acos(vec.cor)*(180/pi)
	return(c(vector.cor=vec.cor, vec.angle=vec.angle))}	

# PLS function from Claude 2010, with modification - updated 4-12-12 WP
PLS <- function( M1, M2 ) {
	p1 <- dim(M1)[2]
	# we allow p2 to be a simple vector, instead of a matrix, thus:-
	if ( is.vector(M1) == T ){
		stop( paste( "'M1' must be specified as a matrix" ) )
		 }
	if ( is.vector(M2) == T ) { p2 <- 1 }
		else { p2 = dim(M2)[2] }
	
	n <- dim(M1)[1]
	sM12 <- svd(var(cbind(M1,M2))[1:p1, (p1+1):(p1+p2)])
	vM12 <- var(cbind(M1,M2))[1:p1, (p1+1):(p1+p2)]
	vM21 <- var(cbind(M1,M2))[(p1+1):(p1+p2), 1:p1]
	v11 <- var(M1)
	v22 <- var(M2)
	D <- sM12$d
	F1 <- sM12$u
	F2 <- sM12$v
	Rv <- sum(diag(vM12%*%vM21))/sqrt(sum(diag(v11%*%v11))*
	sum(diag(v22%*%v22)))
	return( list(Rv_coef=Rv, Singular_Values=D, block1_PLS_vecs=F1, block2_PLS_vecs=F2) )
}


BootstrapRv <- function(data, mat_ind1, mat_ind2, index) {
	data <- data[ index, ] # We will sample along rows of the data frame
	pls_out <- PLS( data[,mat_ind1], data[,mat_ind2] )
	c( pls_out$Rv_coef, pls_out$Singular_Values )	}

RandomizationRv <- function(dat, mat_ind1, mat_ind2){
	dat1 <- dat[ ,mat_ind1]
	dat2 <- dat[sample(nrow(dat), nrow(dat), replace=F), mat_ind2]
	pls_out <- PLS(dat1, dat2)
	c( pls_out$Rv, pls_out$Singular_Values )	}



## ----system_info, echo=FALSE, results='asis'-----------------------------
sessionInfo()


## ----read_in_data_&_data_hygiene, fig.keep='first', echo=FALSE, results='hide'----
######################################
##### read in and munge the data #####
######################################

# NB:– MAKE SURE TO CHANGE THIS LINE TO *YOUR* WORKING DIRECTORY! We currently have it set assuming you start in the directory.
setwd( "./" )

landmarks <- read.csv( "../data/Cricket_Wing_LMs.csv" )
names(landmarks)

calls <- read.csv( "../data/CGCricketCallMeasures.csv" )
names(calls)

length( intersect( landmarks$Male.code, calls$ID ))


# first order of business: flip the left wings so that appear to be right wings
landmarks[ landmarks$Wing=='L', seq( 8, 44, 2 ) ] <- landmarks[ landmarks$Wing=='L', seq( 8, 44, 2 ) ] * -1

dim( landmarks )

# second order of business: fix missing LMs

lmsL <- landmarks[ landmarks$Wing=="L", 7:44]
lmsR <- landmarks[ landmarks$Wing=="R", 7:44]

fixedL <- estimate.missing( arrayspecs( as.matrix( lmsL ) , 19, 2 ), method="TPS" )
fixedR <- estimate.missing( arrayspecs( as.matrix( lmsR ) , 19, 2 ), method="TPS" )

proc_coords <- gpagen( abind(fixedR, fixedL) )$coords
Csize <- gpagen( abind(fixedR, fixedL) )$Csize

ProcCoordsCsize <- data.frame( cbind( two.d.array( proc_coords ), Csize ))
names( ProcCoordsCsize )

ProcCoordsCsize$Digitizing.order <- factor( dimnames( ProcCoordsCsize )[[1]])

landmarks <- merge( landmarks, ProcCoordsCsize, by.x="Digitizing.order", by.y="Digitizing.order" )

names(landmarks)
LM_PCs <- prcomp( landmarks[,45:82] )
summary(LM_PCs)

landmarks <- cbind( landmarks, LM_PCs$x[,1:34] )


######################################
#####    data hygiene etc.		 #######
######################################
plot( landmarks$PC1, landmarks$PC2 )
# identify( landmarks$PC1, landmarks$PC2 )
# 306 307 348 349

landmarks[ c(306, 307, 348, 349) ,]

plot( landmarks[,46], landmarks[,83], type='n', ylim=c(-0.5,0.5), xlim=c(-0.5,0.5) )

for (j in 1:nrow(landmarks)){
for ( i in 1:19 ){
	points( landmarks[ j, seq(45,81,2)[i]], landmarks[ j, seq(46,82,2)[i]] )
}}

j=258
for ( i in 1:19 ){ 	points( landmarks[ j, seq(45,81,2)[i]], landmarks[ j, seq(46,82,2)[i]], col="red" ) }

hist( landmarks$Csize )

# though their Csize seems reasonable, their LM coord.s suggest wings that are not wing-shaped <- clearly there has been 
# some error in digitizing/processing... These four images (representing 2 individuals: LG_AC14 & LG_KL08) will be removed from the dataset
landmarks <- landmarks[ -c(306, 307, 348, 349) ,]


landmarks$side <- factor( substr( landmarks$Male.code, 9, 9 ) )
landmarks$ID <- factor( substr( landmarks$Male.code, 1, 7 ) )

# center centroid size!
mean( landmarks$Csize )		# 20.73964
landmarks$Csize <- scale( landmarks$Csize, scale=F, center=T )



## ----test_for_directional_asymmetry, echo=FALSE, results='hide', fig.keep='first'----
######################################
### test for directional asymmetry ###
######################################


# Directional asymmetry <- does size differ left-to-right?
summary( lm( landmarks$Csize ~ landmarks$side ) )	# no
summary( lm( landmarks$Csize ~ landmarks$side + landmarks$Male.no ) )	# no

plot( landmarks$Csize ~ landmarks$side, xlab="Side", ylab="centroid size (centred)", main="wing size is not asymmetric" )


# Directional asymmetry <- does shape differ left-to-right?
asymm_mod <- lm( as.matrix( landmarks[84:117] ) ~ landmarks$side + landmarks$Csize )
summary( manova( asymm_mod ))	# it would appear yes...
shapePRsq( asymm_mod )	# ...to the tune of 7% of variance

# permutation to make sure...
asymm_mod_perm <- rep( NA, 1000 )
for(i in 1:1000){ 
	asymm_mod_perm[i] <- summary( manova( lm( as.matrix( landmarks[ sample(nrow(landmarks), nrow(landmarks), replace=F) ,84:117] ) ~ landmarks$side + landmarks$Csize ) ))$stats[1,2]}
	hist(asymm_mod_perm, xlim=c(0,1))
	abline( v=summary( manova( asymm_mod ))$stats[1,2], col="red")
	#pseudo-p-val
	mean(c(asymm_mod_perm >= summary( manova( asymm_mod ))$stats[1,2], 1))


# in light of this I'm going to use the mean wing shape within individual
LMs <- aggregate( landmarks[,45:117], by=list( landmarks$ID, landmarks$Male.no, landmarks$Gen., landmarks$Pop.), mean )

names( LMs )[1:4] <- c( "Male.ID", "Male.no", "Gen", "Pop" )

table( LMs[, 3:4] )



## ----testing_for_pop_x_gen_variation_in_size, echo=FALSE, results='hide', fig.keep='last'----
######################################
##### test for variation in size #####
######################################


# SIZE
# does size differ among populations / between rearing env.?
summary( lm( LMs$Csize ~ LMs$Pop * LMs$Gen ) )
summary.aov( lm( LMs$Csize ~ LMs$Pop * LMs$Gen ) )

da.mod1 <- lm( LMs$Csize ~ LMs$Pop * LMs$Gen )
da.mod2 <- lm( LMs$Csize ~ LMs$Pop + LMs$Gen )
anova( da.mod1, da.mod2, test="F" )
# F test indicates that the interaction is important, therefore I ought to use type 3 SS

Anova( lm( LMs$Csize ~ LMs$Pop * LMs$Gen ), type='3', contrasts=list(topic=contr.sum, sys=contr.sum) )

PRsq( lm( LMs$Csize ~ LMs$Pop * LMs$Gen ) )

plot( LMs$Csize ~ LMs$Pop )
plot( LMs$Csize ~ LMs$Gen )
# yes, it appears that size differs <- this fits with the pattern from pronotum width data in the JEB paper...

# pop:gen sizes plot
mod <- lm( Csize ~ Pop * Gen, data= LMs )
# add estimates to the mean value before it was centered
size_est <- summary( Effect( c("Pop", "Gen"), mod, se=T ))$effect + 20.73
size_up <- summary( Effect( c("Pop", "Gen"), mod, se=T ))$upper + 20.73
size_dn <- summary( Effect( c("Pop", "Gen"), mod, se=T ))$lower + 20.73

# pdf( "csize_by_pop.pdf" )
	plot( size_est[,1], ylim=c(18, 25), xlim=c(0.7, 6.7), pch=16, col="blue", xaxt='n', xlab='Population', ylab='forewing centroid size', main="Csize by Pop")
		axis( side=1, seq(1.1, 6.1, 1), c("ACT", "KL", "SA", "SL", "TAS", "WA") )
		points( seq(1.2, 6.2, 1), size_est[,2], pch=16, col="red" )
		legend( 0.7, 25, legend=c( "wild-caught", "lab-reared" ), pch=c(16,16), col=c("blue","red") )

	# plot CI lines
		for (i in 1:6) {
			lines( c( (1:6)[i], (1:6)[i] ), 
					c( size_up[i,1], size_dn[i,1] ), col="blue")
			lines( c( seq(1.2, 6.2, 1)[i], seq(1.2, 6.2, 1)[i] ), 
					c( size_up[i,2], size_dn[i,2] ), col="red")		}
# dev.off()



## ----testing_for_pop_x_gen_variation_in_shape, echo=FALSE, results='hide'----

######################################
##### test for variation in shape ####
######################################

# does shape vary among populations/between rearing env.?
summary( manova( lm( as.matrix( LMs[44:77] ) ~ LMs$Pop * LMs$Gen + LMs$Csize )))

da.mod3 <- lm( as.matrix( LMs[,44:77] ) ~ LMs$Pop + LMs$Gen + LMs$Csize )
da.mod4 <- lm( as.matrix( LMs[,44:77] ) ~ (LMs$Pop + LMs$Gen + LMs$Csize)^2 )
anova( da.mod3, da.mod4 )
# again, the interactions appear important and type 3 SS is appropriate?

shape_mod <- lm( as.matrix( LMs[44:77] ) ~ (LMs$Pop + LMs$Gen + LMs$Csize)^2 )

Manova( shape_mod, test='Wilks', type=3)

# permutation to make sure...
shape_mod_perm <- matrix( NA, nrow=1000, ncol=6 )
for(i in 1:1000){ 
	shape_mod_perm[i,] <- summary( manova( lm( as.matrix( LMs[ sample(nrow(LMs), nrow(LMs), replace=F) ,44:77] ) ~ (LMs$Pop + LMs$Gen + LMs$Csize)^2 ) ))$stats[1:6,2]}
	#pseudo-p-val
	mean(c(shape_mod_perm[,1] >= summary( manova( shape_mod ))$stats[1,2], 1))
	mean(c(shape_mod_perm[,2] >= summary( manova( shape_mod ))$stats[2,2], 1))
	mean(c(shape_mod_perm[,3] >= summary( manova( shape_mod ))$stats[3,2], 1))
	mean(c(shape_mod_perm[,4] >= summary( manova( shape_mod ))$stats[4,2], 1))
	mean(c(shape_mod_perm[,5] >= summary( manova( shape_mod ))$stats[5,2], 1))
	mean(c(shape_mod_perm[,6] >= summary( manova( shape_mod ))$stats[6,2], 1))


shapePRsq( shape_mod )
shapeRsq( lm( as.matrix( LMs[,44:77] ) ~ LMs$Pop ) )
shapeRsq( lm( as.matrix( LMs[,44:77] ) ~ LMs$Gen ) )
shapeRsq( lm( as.matrix( LMs[,44:77] ) ~ LMs$Csize ) )
shapeRsq( lm( as.matrix( LMs[,44:77] ) ~ LMs$Pop:LMs$Gen ) )
shapeRsq( lm( as.matrix( LMs[,44:77] ) ~ LMs$Csize:LMs$Pop ) )
shapeRsq( lm( as.matrix( LMs[,44:77] ) ~ LMs$Csize:LMs$Gen ) )

# let's take the shape residuals to look at further
popgensizeresid <- shape_mod$resid
LMs <- cbind( LMs, popgensizeresid )
names( LMs )[78:111] <- paste( "res.", names( LMs )[78:111], coll='', sep='' )



## ----print_manova_output, echo=FALSE-------------------------------------
Manova( shape_mod, test='Wilks', type=3)


## ----size_&_calls_by_population, echo=FALSE, results='hide', warning=FALSE----
######################################
##### size & calls by population #####
######################################

# at this point we need to combine the call data with the shape data
calls$Male.ID <- factor( substr( calls$ID, 1, 7 ) )
calls <- aggregate( calls[,4:13], list(calls$Male.ID, calls$Generation, calls$Pop.), mean )
names( calls )
names( calls )[1:3] <- c( "Male.ID", "Gen", "Pop" )

call_LM <- merge( LMs, calls, by.x=c("Male.ID"), by.y=c("Male.ID") )

# remove duplicates columns
call_LM <- na.omit( call_LM[, -c(113, 112)] )

# so now, lets PLS size and calls within pop?
call_LM$overall <- factor( call_LM$Pop:call_LM$Gen )
cs_pls_by_pop <- list()
cs_bs_by_pop <- list()
cs_perm_by_pop <- list()

for (i in 1:12) {
	X <- levels( call_LM$overall )[i]
	thiscall <- as.matrix( call_LM[ call_LM$overall==X ,78:111] )
	thissize <- call_LM$Csize[ call_LM$overall==X ]
	cs_pls_by_pop[[i]] <- PLS( thiscall, thissize )
	cs_bs_by_pop[[i]] <- boot( call_LM[ call_LM$overall==X,], BootstrapRv, mat_ind1=78:111, mat_ind2=43, R = 1000 )
	cs_perm_by_pop[[i]] <- replicate(1000, RandomizationRv( call_LM[ call_LM$overall==X,], mat_ind1=78:111, mat_ind2=43 ))
}

# 'p-values' for all axes, all PLS's
for (i in 1:12) {
  print( levels( call_LM$overall )[i] )
	print( mean(c(cs_perm_by_pop[[i]][1,] >= 
					cs_pls_by_pop[[i]][1], 1)) )
	print( mean( c( cs_perm_by_pop[[i]][2,] >= 
					cs_pls_by_pop[[i]][2][[1]], 1)))
}


odds <- seq(1, 11, 2)  ;	evens <- seq(2, 12, 2)

size_pls_vc <- matrix(NA, 12, 12)
colnames(size_pls_vc) <- rownames(size_pls_vc) <- levels(call_LM$overall)
BS_pls_by_pop <- list()
nitt <- 1000

for (i in 1:12) {
	X <- levels( call_LM$overall )[i]
	Xcalls <- call_LM[ call_LM$overall==X ,]
	thiscall <- as.matrix( Xcalls[ ,117:121] )
	thissize <- Xcalls$Csize
	cs_pls_by_pop[[i]] <- PLS( thiscall, thissize )
	BS_pls_by_pop[[i]] <- matrix( NA, ncol=5, nrow=nitt )
	for (k in 1:nitt) {
		BScall <- as.matrix( Xcalls[ sample( 1:dim(Xcalls)[1], dim(Xcalls)[1], rep=T ) ,117:121] )
		BS_pls_by_pop[[i]][k,] <- PLS( BScall, thissize )$block1_PLS_vecs
		}
	}
pls_size_song_by_pop_up <- matrix( NA, ncol= 12, nrow= 12 )
pls_size_song_by_pop_lo <- matrix( NA, ncol= 12, nrow= 12 )

for (i in 1:12) {
	for (j in 1:12) {
		if (i != j) {
			size_pls_vc[i,j] <- ang.vec.abs( cs_pls_by_pop[[i]]$block1_PLS_vecs, cs_pls_by_pop[[j]]$block1_PLS_vecs )[1]
			temp2 <- rep(NA, nitt)
			for (l in 1:nitt) {
				V1b <- BS_pls_by_pop[[i]][l,]
				V2b <- BS_pls_by_pop[[j]][l,]
				temp2[l] <- ang.vec.abs(V1b, V2b)[1]
			}
			pls_size_song_by_pop_up[i,j] <- hpd( temp2 )[1]
			pls_size_song_by_pop_lo[i,j] <- hpd( temp2 )[2]
		}	}	}

size_pls_vc[c(odds),c(evens)]
pls_size_song_by_pop_up[c(odds),c(evens)]
pls_size_song_by_pop_lo[c(odds),c(evens)]



## ----plotting_size_shape_rv, echo=FALSE, results='hide', warning=FALSE----
# now to make a plot of the size of the Rv coefficient among populations

cs_pls_plot_data <- matrix( NA, ncol=5, nrow=12 )
row.names( cs_pls_plot_data ) <- levels( call_LM$overall )

for (i in 1:12) {
	cs_pls_plot_data[i,1] <- mean(cs_bs_by_pop[[i]]$t[,1])
	cs_pls_plot_data[i,2:3] <- boot.ci( cs_bs_by_pop[[i]] , conf = 0.95, type = c("bca", "perc"), index=1)[[4]][1,][4:5]
	cs_pls_plot_data[i,4:5] <- apply( cs_perm_by_pop[[i]], 1, mean )
}

# pdf( "pls_call_size.pdf" )
	plot( cs_pls_plot_data[ seq(1, 11, 2) ,1], ylim=c(0, 1), xlim=c(0.7, 6.7), pch=16, col="blue", xaxt='n', xlab='Population', ylab='RV coefficient', main="PLS with Csize")
	axis( side=1, seq(1.1, 6.1, 1), c("ACT", "KL", "SA", "SL", "TAS", "WA") )
	points( seq(1.2, 6.2, 1), cs_pls_plot_data[ seq(2, 12, 2) ,1], pch=16, col="red" )
	legend( 0.7, 1, legend=c( "wild-caught", "lab-reared" ), pch=c(16,16), col=c("blue","red") )
	# plot CI lines
	odds <- seq(1, 11, 2)	;	evens <- seq(2, 12, 2)
	wild <- 1:6	;	lab <- seq(1.2, 6.2, 1)	; null <- seq(1.3, 6.3, 1)
	for (i in 1:6) {
		lines( c( wild[i], wild[i] ), c( cs_pls_plot_data[odds[i],2], cs_pls_plot_data[odds[i],3] ), col="blue" )
		lines( c( lab[i], lab[i] ), c( cs_pls_plot_data[evens[i],2], cs_pls_plot_data[evens[i],3] ), col="red" )	
		lines( c( null[i], null[i] ), c( cs_pls_plot_data[evens[i],4], cs_pls_plot_data[evens[i],5] ), lwd=4, col="#64646490" )
		}
# dev.off()



## ----shape_&_calls_by_population, echo=FALSE, results='hide', warning=FALSE----
######################################
##### shape & calls by population ####
######################################

pls_by_pop <- list()
bs_by_pop <- list()
perm_by_pop <- list()

	for (i in 1:12) {
		X <- levels( call_LM$overall )[i]
		thiscall <- as.matrix( call_LM[ call_LM$overall==X ,117:121] )
		thisshape <- as.matrix( call_LM[ call_LM$overall==X ,78:111] )
		pls_by_pop[[i]] <- PLS( thisshape, thiscall )
		bs_by_pop[[i]] <- boot( call_LM[ call_LM$overall==X,], BootstrapRv, mat_ind1=117:121, mat_ind2=78:111, R = 1000 )
		perm_by_pop[[i]] <- replicate(1000, RandomizationRv( call_LM[ call_LM$overall==X,], mat_ind1=117:121, mat_ind2=78:111 ))
	}

str(pls_by_pop)
wing_vec_cor_mat <- matrix( NA, ncol= 12, nrow= 12 )
shape_vec_cor_mat <- matrix( NA, ncol= 12, nrow= 12 )

# 'p-values' for all axes, all PLS's
for (i in 1:12) {
  print( levels( call_LM$overall )[i] )
  print( mean(c(perm_by_pop[[i]][1,] >= 
                  pls_by_pop[[i]][1], 1)) )
  print( mean( c( perm_by_pop[[i]][2,] >= 
                    pls_by_pop[[i]][2][[1]], 1)))
}

for (i in 1:12) {
	for (j in 1:12) {
		if (i!=j) {
			wing_i <- pls_by_pop[i][[1]]$block2_PLS_vecs[,1]
			wing_j <- pls_by_pop[j][[1]]$block2_PLS_vecs[,1]
			wing_vec_cor_mat[i,j] <- ang.vec.abs( wing_i, wing_j )[1]
			shape_i <- pls_by_pop[i][[1]]$block1_PLS_vecs[,1]
			shape_j <- pls_by_pop[j][[1]]$block1_PLS_vecs[,1]
			shape_vec_cor_mat[i,j] <- ang.vec.abs( shape_i, shape_j )[1]      
	}	}	}

wing_vec_cor_mat
shape_vec_cor_mat



## ----vector_correlations_between_PLS_vecs_for_pops&gens, echo=FALSE, results='hide', warning=FALSE----
##### bootstrapping the vector correlations between pop.s for the 1st PLS vector from shape vs. calls

nitt <- 1000
pls_shape_song_by_pop_up <- matrix( NA, nrow=12, ncol=12 )
pls_shape_song_by_pop_lo <- matrix( NA, nrow=12, ncol=12 )
all_the_BS <- list()
odds <- seq(1, 11, 2)  ;	evens <- seq(2, 12, 2)

  for (i in 1:12) {
    temp3 <- matrix( NA, nrow=nitt, ncol=34 )
    for (k in 1:nitt) {
      X <- levels( call_LM$overall )[i]
      thiscall <- as.matrix( call_LM[ call_LM$overall==X ,117:121] )
      BS_call <- thiscall[ sample( nrow(thiscall), nrow(thiscall), replace=T ) ,]
      thisshape <- as.matrix( call_LM[ call_LM$overall==X ,78:111] )
      temp3[k,] <- PLS( thisshape, BS_call )$block1[,1]
    }
    all_the_BS[[i]] <- temp3
  }


for (i in 1:12) {
  for (j in 1:12) {
    if (i!=j) {
      temp4 <- rep( NA, nitt )
      for (k in 1:nitt){
        shape_i <- as.matrix( all_the_BS[[i]][k,] )
        shape_j <- as.matrix( all_the_BS[[j]][k,] )
        temp4[k] <- ang.vec.abs( shape_i, shape_j )[1]
        }
      pls_shape_song_by_pop_lo[i,j] <- hpd( temp4 )[1]
      pls_shape_song_by_pop_up[i,j] <- hpd( temp4 )[2]
    }  }	}

wing_vec_cor_mat[c(odds),c(evens)]
pls_shape_song_by_pop_up[c(odds),c(evens)]
pls_shape_song_by_pop_lo[c(odds),c(evens)]



## ----plotting_shape_calls_rv, echo=FALSE, results='hide', fig.keep='last', warning=FALSE----
# now to make a plot of the size of the RV coefficient among populations

pls_plot_data <- matrix( NA, ncol=5, nrow=12 )
row.names( pls_plot_data ) <- levels( call_LM$overall )

for (i in 1:12) {
	pls_plot_data[i,1] <- mean(bs_by_pop[[i]]$t[,1])
	pls_plot_data[i,2:3] <- boot.ci( bs_by_pop[[i]] , conf = 0.95, type = c("bca", "perc"), index=1)[[4]][1,][4:5]
	pls_plot_data[i,4:5] <- hpd( perm_by_pop[[i]][,1] )
}


# pdf( "pls_shape_calls.pdf" )
	plot( pls_plot_data[ seq(1, 11, 2) ,1], ylim=c(0, 1), xlim=c(0.7, 6.7), pch=16, col="blue", xaxt='n', xlab='Population', ylab='RV coefficient', main="PLS with PC residuals")
	axis( side=1, seq(1.1, 6.1, 1), c("ACT", "KL", "SA", "SL", "TAS", "WA") )
	points( seq(1.2, 6.2, 1), pls_plot_data[ seq(2, 12, 2) ,1], pch=16, col="red" )
	legend( 0.7, 1, legend=c( "wild-caught", "lab-reared" ), pch=c(16,16), col=c("blue","red") )
	# plot CI lines
	odds <- seq(1, 11, 2)	;	evens <- seq(2, 12, 2)
	wild <- 1:6	;	lab <- seq(1.2, 6.2, 1)
	for (i in 1:6) {
		lines( c( wild[i], wild[i] ), c( pls_plot_data[odds[i],2], pls_plot_data[odds[i],3] ), col="blue" )
		lines( c( lab[i], lab[i] ), c(pls_plot_data[evens[i],2], pls_plot_data[evens[i],3] ), col="red" )
		lines( c( null[i], null[i] ), c( pls_plot_data[evens[i],4], pls_plot_data[evens[i],5] ), lwd=4, col="#64646490" )	
		}
# dev.off()


## ----printing_vector_correlation_matrix, echo=FALSE----------------------
crap <- wing_vec_cor_mat[c(odds),c(evens)]
row.names( crap ) <- levels( call_LM$overall )[odds]
colnames( crap ) <- levels( call_LM$overall )[evens]
round( crap, 3 )


## ----fossil_code_for_unused_plot, echo=FALSE, results='hide'-------------
##### supplementary plot #####
# pdf( "RV_vs_RV_supp_plot.pdf" )
# plot( pls_plot_data[ seq(1, 11, 2) ,1], cs_pls_plot_data[ seq(1, 11, 2) ,1], pch=16, xlim=c(0,1), ylim=c(0,0.8), col="blue",
#       xlab="RV coef. call vs. shape", ylab="RV coef. call vs size")
#     points( pls_plot_data[ seq(2, 12, 2) ,1], cs_pls_plot_data[ seq(2, 12, 2) ,1], pch=16, col="red" )
#     legend( 0.002, 0.8, legend=c( "wild-caught", "lab-reared" ), pch=c(16,16), col=c("blue","red") )
# 
# for( i in 1:12 ){
#     lines( c( pls_plot_data[odds[i],2], pls_plot_data[odds[i],3] ),
#            c( cs_pls_plot_data[odds[i],1], cs_pls_plot_data[odds[i],1] ), col="blue" )
#     lines( c( pls_plot_data[odds[i],1], pls_plot_data[odds[i],1] ),
#            c( cs_pls_plot_data[odds[i],2], cs_pls_plot_data[odds[i],3] ), col="blue" )
#     lines( c( pls_plot_data[evens[i],2], pls_plot_data[evens[i],3] ),
#            c( cs_pls_plot_data[evens[i],1], cs_pls_plot_data[evens[i],1] ), col="red" )
#     lines( c( pls_plot_data[evens[i],1], pls_plot_data[evens[i],1] ),
#            c( cs_pls_plot_data[evens[i],2], cs_pls_plot_data[evens[i],3] ), col="red" )
#   }
# dev.off()


## ----analysis_of_allometry, echo=FALSE, results='hide', warning=FALSE----
######################################
#####   Allometry by population  #####
######################################

# here I'm going to model the allometry within each Pop/gen and then compare vectors of allometry coefs <- BS FOR CI'S!!!
allom_by_pop_lm <- list()
allom_by_pop_lm_R <- list()
allom_by_pop_Bs <- list()
nitt <- 100
	for (i in  1:12) {
		X <- levels( call_LM$overall )[i]
		Xcalls <- call_LM[call_LM$overall==X, ]
		Xmod <- lm( as.matrix( Xcalls[, 44:77] ) ~ Xcalls$Csize )
		allom_by_pop_lm[[i]] <- coef( Xmod )[2,]
		allom_by_pop_lm_R[[i]] <- shapePRsq( Xmod )
		allom_by_pop_Bs[[i]] <- matrix( NA, ncol=34, nrow=nitt )
		for (j in 1:nitt) {
			BSmod <- lm( as.matrix( Xcalls[ sample( 1:dim(Xcalls)[1], dim(Xcalls)[1], rep=T ) , 44:77] ) ~ Xcalls$Csize )
			allom_by_pop_Bs[[i]][j,] <- coef( BSmod )[2,]
		}
	}
allom_by_pop_lm_VC <- matrix( NA, ncol= 12, nrow= 12 )
allom_by_pop_VC_up <- matrix( NA, ncol= 12, nrow= 12 )
allom_by_pop_VC_lo <- matrix( NA, ncol= 12, nrow= 12 )
	for (j in 1:12) {
		for (k in 1:12) {
			if (j!=k) {
				V1 <- allom_by_pop_lm[[j]]
				V2 <- allom_by_pop_lm[[k]]
				allom_by_pop_lm_VC[j,k] <- ang.vec.abs(V1, V2)[1]
				temp <- rep( NA, nitt )
				for (l in 1:nitt) {
					V1b <- allom_by_pop_Bs[[j]][l,]
					V2b <- allom_by_pop_Bs[[k]][l,]
					temp[l] <- ang.vec.abs(V1b, V2b)[1]
				}
				allom_by_pop_VC_up[j,k] <- hpd( temp )[1]
				allom_by_pop_VC_lo[j,k] <- hpd( temp )[2]
			}	}	}

colnames( allom_by_pop_lm_VC ) <- levels( call_LM$overall )
rownames( allom_by_pop_lm_VC ) <- levels( call_LM$overall )

allom_by_pop_lm_VC[c(odds),c(evens)]
allom_by_pop_VC_up[c(odds),c(evens)]
allom_by_pop_VC_lo[c(odds),c(evens)]


## ----print_allometry_VC_matrix, echo=FALSE-------------------------------
round( allom_by_pop_lm_VC[c(odds),c(evens)], 3 )


## ----discriminant_function_plot, fig.keep='last', echo=FALSE, results='hide'----
######################################
#####   LDA by population/gen.   #####
######################################

pop_gen <- lda( call_LM[,44:77], call_LM$Pop.x:call_LM$Gen.x )

pop_gen_scores <- as.matrix(call_LM[,44:77]) %*% as.matrix(pop_gen$scaling)

call_LM <- data.frame( call_LM, pop_gen_scores )


df_plot_data <- ddply( call_LM, .(Pop.x, Gen.x), summarize, mLD1=mean( LD1 ), mLD2=mean( LD2 ), sdLD1=sd( LD1 ), sdLD2=sd( LD2 ) )

df_plot_data$pch <- rep( c(22, 24), 6 )
popcols <- c("red", "blue", "yellow", "black", "orange", "green")
df_plot_data$col <- c( interleave( popcols, popcols ))

# pdf( "wing_form_DF_plot.pdf" )
	plot( call_LM$LD1, call_LM$LD2, xlab='1st linear discriminant', ylab='2nd linear discriminant', col="#64646450", pch=16, xlim=c(-6,6) )	
	for (i in 1:nrow(df_plot_data)) {
		lines( c( df_plot_data$mLD1[i] + df_plot_data$sdLD1[i],
				  df_plot_data$mLD1[i] - df_plot_data$sdLD1[i] ),
				c( df_plot_data$mLD2[i], df_plot_data$mLD2[i] ) )
		lines( c( df_plot_data$mLD1[i], df_plot_data$mLD1[i] ),
				c( df_plot_data$mLD2[i] + df_plot_data$sdLD2[i],
			   	  df_plot_data$mLD2[i] - df_plot_data$sdLD2[i] ) )
		}
	points( df_plot_data[,3:4], pch=df_plot_data$pch, bg=df_plot_data$col, cex=1.5 )
# dev.off()


## ----between_group_covariance, echo=FALSE, results='hide', fig.keep='none'----

# Compute population means for calls & shape – look at RV coef.
popmeansF <- aggregate( call_LM[ call_LM$Gen.x=="FIELD" , c( 44:77, 117:121 )], by= list( "Pop"=call_LM$Pop.x[ call_LM$Gen.x=="FIELD"] ) , mean )  #check call cols

popmeansL <- aggregate( call_LM[ call_LM$Gen.x=="LAB" , c( 44:77, 117:121 )], by= list( "Pop"=call_LM$Pop.x[ call_LM$Gen.x=="LAB"] ) , mean )  #check call cols


pop_mean_plsF <- PLS( as.matrix( popmeansF[, 2:35] ), as.matrix( popmeansF[, 36:40] ) )
pop_mean_plsL <- PLS( as.matrix( popmeansL[, 2:35] ), as.matrix( popmeansL[, 36:40] ) )

# permute all possible combinations (720 = 6 factorial) to compare with this RV
all_perms <- permn( 1:6 )
pop_permsF <- rep( NA, 720 )
pop_permsL <- rep( NA, 720 )

for( i in 1:factorial( 6 )) {
      thisperm <- c(1:6)[ all_perms[[i]] ]
      thisplsF <- PLS( as.matrix( popmeansF[, 2:35] ), as.matrix( popmeansF[ thisperm, 36:40] ) )
      thisplsL <- PLS( as.matrix( popmeansL[, 2:35] ), as.matrix( popmeansL[ thisperm, 36:40] ) )
      pop_permsF[i] <- thisplsF$Rv_coef
      pop_permsL[i] <- thisplsL$Rv_coef
  }

hist( c(pop_permsF, pop_permsL), xlim=c(0, 1) )   ;   abline( v=pop_mean_plsF$Rv_coef, col="red", lwd=2 )
hist( c(pop_permsF, pop_permsL), xlim=c(0, 1) )   ;   abline( v=pop_mean_plsL$Rv_coef, col="red", lwd=2 )

length( pop_permsF[ pop_permsF > pop_mean_plsF$Rv_coef ] ) / factorial( 6 )
length( pop_permsL[ pop_permsL > pop_mean_plsL$Rv_coef ] ) / factorial( 6 )



## ----comparing_the_major_PLS_axes_between_pops, results='hide', echo=FALSE, fig.keep='none'----

betweenF1 <- pop_mean_plsF$block1_PLS_vecs
betweenF2 <- pop_mean_plsF$block2_PLS_vecs

betweenL1 <- pop_mean_plsL$block1_PLS_vecs
betweenL2 <- pop_mean_plsL$block2_PLS_vecs

within1 <- list( pls_by_pop[[1]][3], pls_by_pop[[2]][3], pls_by_pop[[3]][3], pls_by_pop[[4]][3], 
                pls_by_pop[[5]][3], pls_by_pop[[6]][3], pls_by_pop[[7]][3], pls_by_pop[[8]][3],
                pls_by_pop[[9]][3], pls_by_pop[[10]][3], pls_by_pop[[11]][3], pls_by_pop[[12]][3] )
within2 <- list( pls_by_pop[[1]][4], pls_by_pop[[2]][4], pls_by_pop[[3]][4], pls_by_pop[[4]][4], 
                pls_by_pop[[5]][4], pls_by_pop[[6]][4], pls_by_pop[[7]][4], pls_by_pop[[8]][4],
                pls_by_pop[[9]][4], pls_by_pop[[10]][4], pls_by_pop[[11]][4], pls_by_pop[[12]][4] )

within_withoutF1 <- rep(NA, 12)  ;  within_withoutF2 <- rep(NA, 12)
within_withoutL1 <- rep(NA, 12)  ;  within_withoutL2 <- rep(NA, 12)

for (i in 1:12) {
  within_withoutF1[i]  <- ang.vec.abs( betweenF1[,1], within1[[i]]$block1_PLS_vecs[,1] )[1]
  within_withoutF2[i] <- ang.vec.abs( betweenF2[,1], within2[[i]]$block2_PLS_vecs[,1] )[1]
  within_withoutL1[i]  <- ang.vec.abs( betweenL1[,1], within1[[i]]$block1_PLS_vecs[,1] )[1]
  within_withoutL2[i] <- ang.vec.abs( betweenL2[,1], within2[[i]]$block2_PLS_vecs[,1] )[1]
}

summary( within_withoutF1[ seq(1, 11, 2) ] )
summary( within_withoutL1[ seq(2, 12, 2) ] )

summary( within_withoutF2[ seq(1, 11, 2) ] )
summary( within_withoutL2[ seq(2, 12, 2) ] )



## ----wing_shape_variation_plot, results='hide', echo=FALSE, fig.height=5, fig.width=5----

meanshape <- colMeans( call_LM[, 5:42] )

AC_F <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="FIELD" & call_LM$Pop.x=="AC", 5:42] )
KL_F <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="FIELD" & call_LM$Pop.x=="KL", 5:42] )
SA_F <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="FIELD" & call_LM$Pop.x=="SA", 5:42] )
SL_F <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="FIELD" & call_LM$Pop.x=="SL", 5:42] )
TA_F <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="FIELD" & call_LM$Pop.x=="TA", 5:42] )
WA_F <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="FIELD" & call_LM$Pop.x=="WA", 5:42] )

AC_L <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="LAB" & call_LM$Pop.x=="AC", 5:42] )
KL_L <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="LAB" & call_LM$Pop.x=="KL", 5:42] )
SA_L <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="LAB" & call_LM$Pop.x=="SA", 5:42] )
SL_L <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="LAB" & call_LM$Pop.x=="SL", 5:42] )
TA_L <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="LAB" & call_LM$Pop.x=="TA", 5:42] )
WA_L <- meanshape - colMeans( call_LM[ call_LM$Gen.x=="LAB" & call_LM$Pop.x=="WA", 5:42] )

linkmat <- matrix( c(1,2,2,3,3,4,1,19,19,18,18,17,17,16,11,9,9,8,8,7,7,6,6,5,11,10,10,13,13,12,13,15,15,12,12,11
                     ,12,15,10,6,14,4,14,16,11,17,8,2,6,4,3,6), ncol=2, byrow = T )

popgen <- list( AC_F, KL_F, SA_F, SL_F, TA_F, WA_F, AC_L, KL_L, SA_L, SL_L, TA_L, WA_L )
popgencol <- c( rep( "grey", 6 ), rep( "black", 6 ) )
odds <- seq(1, 37, 2)  ;  evens <- seq(2, 38, 2)  ;   scale.fac <- 5


popgencol <- c( rep( c( "red", "blue", "yellow", "black", "orange", "dark green" ), 2 ) )

# pdf( "pretty_wing_pic_F.pdf", width = 5, height = 5 )
  plot( -0.5:1.2, -0.5:1.2, type='n', xlab='', ylab='', xaxt='n', yaxt='n', main="wing shape mean - FIELD" )
    for (i in 1:6) {
        thispop <- meanshape + ( unlist(popgen[i]) * scale.fac )
        LMx <- thispop[ c(odds) ]
        LMy <- thispop[ c(evens) ]
        thiscol <- popgencol[i]
  for (k in 1:26) {
          lines( c( LMx[ linkmat[k,1] ], LMx[ linkmat[k,2] ] ), 
                 c( LMy[ linkmat[k,1] ], LMy[ linkmat[k,2] ] ), 
                 lwd=0.75, col="dark grey" )
        }
  for (j in 1:19) {
            points( LMx[j], LMy[j], pch=21, cex=0.75, bg=thiscol, col="dark grey" )
        }
      }
# dev.off()

# pdf( "pretty_wing_pic_L.pdf", width = 5, height = 5 )
  plot( -0.5:1.2, -0.5:1.2, type='n', xlab='', ylab='', xaxt='n', yaxt='n', main="wing shape mean - LAB" )
    for (i in 7:12) {
        thispop <- meanshape + ( unlist(popgen[i]) * scale.fac)
        LMx <- thispop[ c(odds) ]
        LMy <- thispop[ c(evens) ]
        thiscol <- popgencol[i]
  for (k in 1:26) {
          lines( c( LMx[ linkmat[k,1] ], LMx[ linkmat[k,2] ] ), 
                 c( LMy[ linkmat[k,1] ], LMy[ linkmat[k,2] ] ), 
                 lwd=0.75, col="dark grey" )
        }
  for (j in 1:19) {
            points( LMx[j], LMy[j], pch=21, cex=0.75, bg=thiscol, col="dark grey" )
        }
      }
# dev.off()



