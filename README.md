# UdacityRadar
Udacity Nanodegree Self Driving Car/ Radar proj


# Project Overview
![ProjectOverview](https://github.com/Felix-yuan2018/SFND_Radar_Target_Generation_and_Detection/blob/master/img/ProjectOverview.png)

* Configure the FMCW waveform based on the system requirements.
* Define the range and velocity of target and simulate its displacement.
* For the same simulation loop process the transmit and receive signal to determine the beat signal
* Perform Range FFT on the received signal to determine the Range
* Towards the end, perform the CFAR processing on the output of 2nd FFT to display the target.

## 1. Radar Specifications
* Frequency of operation = 77 GHz
* Max Range = 200 m
* Range Resolution = 1 m
* Max Velocity = 70 m/s

```matlab
fc= 77e9;             %carrier freq
c=3e8;                %speed of light
max_range=200;        %max range
range_resolution=1;   % range resolution can determine bandwith
max_velocity=100;     % max velocity can determine Ts_x-axis

```


## 2. user Defined Range and Velocity of target

* define the target's initial posotion and velocity
* Note: Velocity remains constant

```matlab
Range = 110;     % Initial position of the target
velocity = -20;  % Initial Velocity of the target

```


## 3. FMCW Waveform Generation

* Design the FMCW waveform by giving the specs of each of its parameters.
* Calculate the Bandwidth (Bsweep), Chirp Time (Tchirp) and Slope (slope) of the FMCW chirp using the requirements above.
* Operating carrier frequency of Radar

```matlab
Bandwidth= c/(2*range_resolution);  
Tchirp = 5.5*2*max_range/c;      %Chirp Time
slope=Bandwidth/Tchirp;          %slope of the FMCW
```

* The number of chirps in one sequence.
* Its ideal to have 2^value for the ease of running the FFT for Doppler Estimation.

```matlab
Nd=128; 
```

* The number of samples on each chirp.

```matlab
Nr=1024;
```

* Timestamp for running the displacement scenario for every sample on each chirp

```matlab
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
```

* Creating the vectors for Tx, Rx and mix based on the total samples input.

```matlab
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal
```

* Similar vectors for range_covered and time delay.

```matlab
r_t=zeros(1,length(t));
td=zeros(1,length(t));
```


## 4. Signal generation and Moving Target simulation


```matlab
for i=1:length(t)    
    
    % For each time stamp update the Range of the Target for constant velocity. 
     r_t(i) = range + velocity*t(i);      %Update the Range of the Target for constant velocity.
     td(i) = (2*r_t(i))/c;                  %compute time_delay.


    
    % For each time sample we need update the transmitted and Received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+ 0.5*slope*t(i)^2));
    Rx(i)  = cos(2*pi*(fc*(t(i)-td(i))+ 0.5*slope*(t(i)-td(i))^2));

    
    % Now by mixing the Transmit and Receive generate the beat signal This is done by element
    % wise matrix multiplication of Transmit and Receiver Signal
    Mix(i) = Tx(i).*Rx(i);
    
end
```


## 5. Range Measurement

* Reshape the vector into Nr*Nd array.
* Nr and Nd here would also define the size of Range and Doppler FFT respectively.
* run the FFT on the beat signal along the range bins dimension (Nr) and normalize
* Take the absolute value of FFT output

```matlab
Mix = reshape(Mix, [Nr,  Nd]);
%run the FFT on the beat signal along the range bins dimension (Nr) and
sig_fft1 = fft(Mix,Nr); 
 
%normalize.
sig_fft1 = sig_fft1./Nr;

% Take the absolute value of FFT output
sig_fft1 = abs(sig_fft1);  


```

* Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
* Hence we throw out half of the samples.

```matlab
single_side_sig_fft1 = sig_fft1(1:Nr/2);

```

* Plotting the range, plot FFT output

```matlab
figure ('Name','Range from First FFT')
plot(single_side_sig_fft1); 
axis ([0 200 0 1]);
```

* Simulation Result

![RangeFFT](https://github.com/Felix-yuan2018/SFND_Radar_Target_Generation_and_Detection/blob/master/img/RangeFFT.png)


## 6. Range Doppler Response

* The 2D FFT implementation is already provided here.
* This will run a 2DFFT on the mixed signal (beat signal) output and generate a range doppler map.
* You will implement CFAR on the generated RDM Range Doppler Map Generation.
* The output of the 2D FFT is an image that has reponse in the range and doppler FFT bins.
* So, it is important to convert the axis from bin sizes to range and doppler based on their Max values.

```matlab
Mix=reshape(Mix,[Nr,Nd]);
```

* 2D FFT using the FFT size for both dimensions.

```matlab
sig_fft2 = fft2(Mix,Nr,Nd);
```

* Taking just one side of signal from Range dimension.

```matlab
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;
```

* Use the surf function to plot the output of 2DFFT and to show axis in both dimensions

```matlab
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
```

* Simulation Result

![DopplerFFT](https://github.com/Felix-yuan2018/SFND_Radar_Target_Generation_and_Detection/blob/master/img/DopplerFFT.png)


## 7. CFAR implementation

* Slide Window through the complete Range Doppler Map
* Select the number of Training Cells in both the dimensions.

```matlab
Td = 10; 
Tr = 8;
```

* Select the number of Guard Cells in both dimensions around the Cell under test (CUT) for accurate estimation

```matlab
Gd = 4; 
Gr = 4;
```

* Offset the threshold by SNR value in dB

```matlab
offset = 6;
```

* Create a vector to store noise_level for each iteration on training cells

```matlab
noise_level = zeros(1,1);
```

* Create a vector to store noise_level for each iteration on training cells design a loop such that it slides the CUT across range doppler map by giving margins at the edges for Training and Guard Cells.
* For every iteration sum the signal level within all the training cells.
* To sum convert the value from logarithmic to linear using db2pow function.
* Average the summed values for all of the training cells used.
* After averaging convert it back to logarithimic using pow2db.
* Further add the offset to it to determine the threshold.
* Next, compare the signal under CUT with this threshold.
* If the CUT level > threshold assign % it a value of 1, else equate it to 0.
* Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR

```matlab
RESULT=zeros(Nr/2,Nd);
 
for i= Tr+Gr+1: Nr/2-(Gr+Tr)
    for j= Td+Gd+1:Nd-(Gd+Td)
        noise_level=0;
        train_cell=0;
        %sliding the grid
        for p= i-(Tr+Gr):i+Tr+Gr
            for q= j-(Td+Gd):j+Td+Gd
                %Verify that there is Guard Cells in both dimensions around the Cell under test (CUT) 
                if(abs(i-p)>Gr||abs(j-q)>Gd)
                    noise_level=noise_level+db2pow(RDM(p,q));
                    train_cell=train_cell+1;
                end
            end
        end
        
        noise_average = noise_level/train_cell;
        noise_threshold = pow2db(noise_average);
        noise_threshold = noise_threshold+offset;
        
        CUT=RDM(i,j);
        if(CUT<noise_threshold)
            RESULT(i,j)=0;
        else
            RESULT(i,j)=1;
        end
    end
end
```

* Display the CFAR output using the Surf function like we did for Range
* Doppler Response output.

```matlab
figure,surf(doppler_axis,range_axis,RESULT);
colorbar;
```

* Simulation Result

![CA-CFAR_Filtered_RDM](https://github.com/Felix-yuan2018/SFND_Radar_Target_Generation_and_Detection/blob/master/img/CA-CFAR_Filtered_RDM.png)
