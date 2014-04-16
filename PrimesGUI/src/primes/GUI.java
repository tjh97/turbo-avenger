package primes;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import javax.sound.sampled.*;
import javax.swing.*;
import javax.swing.event.*;

public class GUI extends JFrame implements Runnable {
	private JTextField textField;

	
	public GUI() {
//		JFrame mainWindow= new JFrame("Prime Calculator");
//		
//		initiateGUI(mainWindow);
//		
		
	}
	
	private void initiateGUI(JFrame mainWindow) {
		mainWindow.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		JPanel panel = new JPanel();
		mainWindow.getContentPane().add(panel, BorderLayout.NORTH);
		
		textField = new JTextField();
		panel.add(textField);
		textField.setColumns(10);
		
		JButton btnIsItPrime = new JButton("Is it Prime?");
		panel.add(btnIsItPrime);
		
		
		pack();
		
		mainWindow.setVisible(true);
		panel.setVisible(true);
		
	}
	
	public void run() {
		
		while (true) {
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
	        System.out.print("Enter Integer:");
	        try{
	            long i = Long.parseLong(br.readLine());
	            if (Prime.isPrime(i)) {
	            	System.out.println("Yes!");
	            } else {
	            	System.out.println("No!");
	            }
	        }catch(NumberFormatException nfe){
	            System.err.println("Invalid Format!");
	        }catch(IOException iofe) {
	        	System.err.println("Invalid Format!");
	        }
	        
		}
	}

}
