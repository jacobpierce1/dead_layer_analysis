import libjacob.meas as meas 


x = meas.meas( 1, 1 )
# print x + x 
# print x * 2
# print 2 * x 
# print 1 / x
# print 1 / meas.cos( x )
y =  meas.meas.from_list( [ meas.meas( 1, 2), meas.meas(3,3) ] )
print y[1] 
