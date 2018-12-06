I'm looking at the k24 vs gamma space, similar to in the polymorphism paper, except with omega = Lambda != 0.

Yesterday (2018-12-05), I attempted this calculation, but ran out of time for omega=Lambda=1 and omega=Lambda=3. My code broke for omega=Lambda=0 because I forgot to reset the next guess of eta and delta if they are NAN. My code broke for omega=Lambda=10.0 and omega=Lambda=30.0 because the radius was being driven above 20 for large gamma values (as expected), which is where I have my sharp cutoff limits (I set E=1e100 if R<0 or R>20).

Today, I will do similar calculations with only two values of k24, 0 and 0.01, to see how long it takes approximately to raster across from gamma=0.01 to gamma=0.2 with 96 steps. I will increase the time allowed for computation to fix the issue mentioned above for the omega=Lambda=1 and omega=Lambda=3 cases.

I will deal with the omega=Lambda=0 case (by resetting eta and delta correctly) and the large omega and Lambda cases (by just quitting the current calculation vs exiting to system) soon.