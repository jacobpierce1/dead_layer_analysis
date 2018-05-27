import matplotlib 
import time

matplotlib.use('Agg')
import matplotlib.pyplot as plt 

for i in range( 10 ) :

    told = time.time() 

    # time.sleep(3)

    f, axarr = plt.subplots( 2, 2 , figsize = (20,10) )

    axarr[0][1].plot( [1,2,3], [4,5,6] )

    plt.savefig( 'test.png', format = 'png' ) 

    tnew = time.time()

    print( tnew - told)
    told = tnew 
    
    plt.cla() 
    plt.close( f ) 
