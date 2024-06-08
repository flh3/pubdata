
data(pisa2012, package = 'MLMusingR') 

system.time({
  m1 <- mixPV(pv1math + pv2math + pv3math + pv4math + pv5math ~
                st29q03 + sc14q02 + st04q01 + escs + (escs|schoolid), 
              weights = c('w_fstuwt', 'w_fschwt'), 
              data = pisa2012)
})



system.time({
  m1 <- mixPV(pv1math + pv2math + pv3math + pv4math + pv5math ~
                st29q03 + sc14q02 + st04q01 + escs + (escs|schoolid), 
              weights = c('w_fstuwt', 'w_fschwt'), 
              data = pisa2012, mc = TRUE)
})

st <- Sys.time()
m1a <- mixPV(pv1math + pv2math + pv3math + pv4math + pv5math ~
              st29q03 + sc14q02 + st04q01 + escs + (1|schoolid), 
            weights = c('w_fstuwt', 'w_fschwt'), 
            data = pisa2012, mc = FALSE)
Sys.time() - st

st <- Sys.time()
m1 <- mixPV(pv1math + pv2math + pv3math + pv4math + pv5math ~
              st29q03 + sc14q02 + st04q01 + escs + (1|schoolid), 
            weights = c('w_fstuwt', 'w_fschwt'), 
            data = pisa2012, mc = TRUE)
Sys.time() - st

summary(m1a)
summary(m1)
