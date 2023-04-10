using System;

namespace ChainExample
{
    class Program
    {
        static void Main(string[] args)
        {
            int numSamples = 10000; // Number of samples to generate
            double mean = 0.0; // Mean of the Gaussian distribution
            double stdDev = 1.0; // Standard deviation of the Gaussian distribution

            // Initial state of the Markov chain
            double currentState = 0.0;

            // Create a random number generator
            Random random = new Random();

            // Array to store the samples
            double[] samples = new double[numSamples];

            // Generate samples using Metropolis-Hastings algorithm
            for (int i = 0; i < numSamples; i++)
            {
                // Propose a new state by adding a random value from a Gaussian distribution
                double proposedState = currentState + RandomGaussian(random) * stdDev;

                // Calculate the acceptance ratio
                double acceptanceRatio = GaussianDensity(proposedState, mean, stdDev) / GaussianDensity(currentState, mean, stdDev);

                // Generate a uniform random number
                double u = random.NextDouble();

                // Accept or reject the proposed state based on the acceptance ratio
                if (u < acceptanceRatio)
                {
                    currentState = proposedState;
                }

                // Store the current state as a sample
                samples[i] = currentState;
            }

            // Output the samples
            Console.WriteLine("Samples from the Gaussian distribution:");
            for (int i = 0; i < numSamples; i++)
            {
                Console.WriteLine(samples[i]);
            }
        }

        // Helper function to generate a random number from a Gaussian distribution
        static double RandomGaussian(Random random)
        {
            double u1 = 1.0 - random.NextDouble();
            double u2 = 1.0 - random.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Sin(2.0 * Math.PI * u2);
            return randStdNormal;
        }

        // Helper function to calculate the density of a Gaussian distribution
        static double GaussianDensity(double x, double mean, double stdDev)
        {
            double exponent = -0.5 * ((x - mean) / stdDev) * ((x - mean) / stdDev);
            return Math.Exp(exponent) / (stdDev * Math.Sqrt(2.0 * Math.PI));
        }
    }
}
