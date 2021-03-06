# function for returning instring without extension (i.e. ".txt" or ".trees")
rmext = function (instring, delimiter="\\.")
{
    library('stringr')
    # get all instances of "."
    alist = str_locate_all(instring, delimiter)
    # get position of last element
    pos2 = tail(alist[[1]],1)[1]-1
    return(substr(instring, 1, pos2))
}

# function for returning a list of letters from a string
str2lst = function (instring)
{
    unlist(strsplit(instring,""))
}

# function to check if return is integer(0)
is.integer0 = function (instring)
{
  is.integer(instring) && length(instring) == 0L
}

# function to check if return is character(0)
is.character0 = function (instring)
{
  is.character(instring) && length(instring) == 0
}

# function for returning listcolumn x in a matrix of lists ("ReTurN List Column")
rtnlc = function (m,x)
{
    aFun = function(r,c,m,x) {as.numeric(unlist(strsplit(m[r,c],",")))[x]}
    vect_aFun = Vectorize(aFun,vectorize.args = c('r','c'))
    return(outer(1:nrow(m), 1:ncol(m) , FUN = vect_aFun, m, x))
}

