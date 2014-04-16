package primes;

import java.util.*;

public class Prime implements Runnable {
	
	public static ArrayList<Long> primeList= new ArrayList<Long>();
	public static Long largestPrime= new Long(2);
	
	
	/** Constructor: num is a prime number.
	 * 
	 * @param num num is a prime number that is the next prime number
	 * in the current list of primes.
	 */
	Prime () {
	}
	
	public static boolean isPrime(long num) {
		synchronized (primeList) {
			if (primeList.contains(num)) return true;
		}
		
		return calculateIsPrime(num);
	}
	
	private static boolean calculateIsPrime(long num) {
		long upperLimit= (long) Math.sqrt(num) + 1;	//last number that needs to be checked
		long prime;	//store prime number from the array
		
		for (int i= 0; i < primeList.size()-1; i++) {
			synchronized (primeList) {
				prime= primeList.get(i);
			}
			if (upperLimit < prime) return true;
			if (num%prime == 0) return false;
		}
		
		return true;
	}
	
	public static void addPrime(long num) {
		// Add new prime to the prime list
		synchronized (primeList) {
			primeList.add(num);
		}
		// Update the value of the largest prime
		synchronized (largestPrime) {
			largestPrime= num;
		}
		
	}
	
	static void calculatePrimes() {
		primeList.add(new Long(2));
		long number= 3;
		
		while (true) {
			long upperLimit= (long) Math.sqrt(number) + 1;	//last number that needs to be checked
			long prime;	//store prime number from the array
			
			for (int i= 0; i < primeList.size(); i++) {
				synchronized (primeList) {
					prime= primeList.get(i);
				}
				
				if (number%prime == 0) break;
				if (upperLimit <= prime) {
					addPrime(number);
					//System.out.println(number);
					break;
				}
			}
			
			number += 2;
		}
		
	}
	
	@Override
	public void run () {
		calculatePrimes();
	}
}
