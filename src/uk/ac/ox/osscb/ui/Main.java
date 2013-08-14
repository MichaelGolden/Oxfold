package uk.ac.ox.osscb.ui;

import java.awt.BorderLayout;
import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.BoxLayout;
import javax.swing.border.TitledBorder;
import javax.swing.JButton;
import javax.swing.SwingConstants;
import java.awt.FlowLayout;
import javax.swing.JTextField;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.Insets;
import javax.swing.JLabel;
import javax.swing.JPasswordField;
import javax.swing.JRadioButton;
import com.jgoodies.forms.layout.FormLayout;
import com.jgoodies.forms.layout.ColumnSpec;
import com.jgoodies.forms.layout.RowSpec;
import com.jgoodies.forms.factories.FormFactory;
import javax.swing.UIManager;
import javax.swing.JSpinner;
import java.awt.GridLayout;
import javax.swing.SpinnerNumberModel;

public class Main extends JFrame {

	private JPanel contentPane;
	private JTextField textField;
	private JTextField textField_1;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					Main frame = new Main();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public Main() {
		setTitle("Oxfold");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 501, 407);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(new BoxLayout(contentPane, BoxLayout.Y_AXIS));
		
		JPanel panel_6 = new JPanel();
		contentPane.add(panel_6);
		panel_6.setLayout(new BoxLayout(panel_6, BoxLayout.Y_AXIS));
		
		JPanel panel = new JPanel();
		panel_6.add(panel);
		panel.setBorder(new TitledBorder(null, "Select an alignment file", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		panel.setLayout(new BorderLayout(0, 0));
		
		JPanel panel_3 = new JPanel();
		panel.add(panel_3, BorderLayout.CENTER);
		GridBagLayout gbl_panel_3 = new GridBagLayout();
		gbl_panel_3.columnWidths = new int[]{0, 0, 0};
		gbl_panel_3.rowHeights = new int[]{0, 0, 0, 0, 0, 0};
		gbl_panel_3.columnWeights = new double[]{0.0, 1.0, Double.MIN_VALUE};
		gbl_panel_3.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
		panel_3.setLayout(gbl_panel_3);
		
		JButton btnNewButton = new JButton("Browse..");
		GridBagConstraints gbc_btnNewButton = new GridBagConstraints();
		gbc_btnNewButton.insets = new Insets(0, 0, 5, 5);
		gbc_btnNewButton.gridx = 0;
		gbc_btnNewButton.gridy = 0;
		panel_3.add(btnNewButton, gbc_btnNewButton);
		
		textField = new JTextField();
		GridBagConstraints gbc_textField = new GridBagConstraints();
		gbc_textField.insets = new Insets(0, 0, 5, 0);
		gbc_textField.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField.gridx = 1;
		gbc_textField.gridy = 0;
		panel_3.add(textField, gbc_textField);
		textField.setColumns(10);
		
		JLabel lblFormat = new JLabel("File format");
		GridBagConstraints gbc_lblFormat = new GridBagConstraints();
		gbc_lblFormat.anchor = GridBagConstraints.EAST;
		gbc_lblFormat.insets = new Insets(0, 0, 5, 5);
		gbc_lblFormat.gridx = 0;
		gbc_lblFormat.gridy = 1;
		panel_3.add(lblFormat, gbc_lblFormat);
		
		JLabel lblNewLabel_2 = new JLabel("Various alignment formats are recognized (e.g. fasta, phylip, nexus)");
		GridBagConstraints gbc_lblNewLabel_2 = new GridBagConstraints();
		gbc_lblNewLabel_2.anchor = GridBagConstraints.WEST;
		gbc_lblNewLabel_2.insets = new Insets(0, 0, 5, 0);
		gbc_lblNewLabel_2.gridx = 1;
		gbc_lblNewLabel_2.gridy = 1;
		panel_3.add(lblNewLabel_2, gbc_lblNewLabel_2);
		
		JLabel lblNewLabel = new JLabel("# sequences");
		lblNewLabel.setHorizontalAlignment(SwingConstants.LEFT);
		GridBagConstraints gbc_lblNewLabel = new GridBagConstraints();
		gbc_lblNewLabel.anchor = GridBagConstraints.EAST;
		gbc_lblNewLabel.insets = new Insets(0, 0, 5, 5);
		gbc_lblNewLabel.gridx = 0;
		gbc_lblNewLabel.gridy = 2;
		panel_3.add(lblNewLabel, gbc_lblNewLabel);
		
		JLabel lblNewLabel_1 = new JLabel("Avg. length");
		GridBagConstraints gbc_lblNewLabel_1 = new GridBagConstraints();
		gbc_lblNewLabel_1.anchor = GridBagConstraints.EAST;
		gbc_lblNewLabel_1.insets = new Insets(0, 0, 5, 5);
		gbc_lblNewLabel_1.gridx = 0;
		gbc_lblNewLabel_1.gridy = 3;
		panel_3.add(lblNewLabel_1, gbc_lblNewLabel_1);
		
		JLabel lblLength = new JLabel("Total length");
		GridBagConstraints gbc_lblLength = new GridBagConstraints();
		gbc_lblLength.anchor = GridBagConstraints.EAST;
		gbc_lblLength.insets = new Insets(0, 0, 0, 5);
		gbc_lblLength.gridx = 0;
		gbc_lblLength.gridy = 4;
		panel_3.add(lblLength, gbc_lblLength);
		
		JPanel panel_1 = new JPanel();
		panel_6.add(panel_1);
		panel_1.setBorder(new TitledBorder(UIManager.getBorder("TitledBorder.border"), "Select a phylogenetic tree", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		panel_1.setLayout(new BorderLayout(0, 0));
		
		JPanel panel_4 = new JPanel();
		panel_1.add(panel_4, BorderLayout.NORTH);
		panel_4.setLayout(new FormLayout(new ColumnSpec[] {
				FormFactory.RELATED_GAP_COLSPEC,
				FormFactory.DEFAULT_COLSPEC,
				FormFactory.RELATED_GAP_COLSPEC,
				ColumnSpec.decode("default:grow"),},
			new RowSpec[] {
				FormFactory.RELATED_GAP_ROWSPEC,
				FormFactory.DEFAULT_ROWSPEC,
				FormFactory.RELATED_GAP_ROWSPEC,
				FormFactory.DEFAULT_ROWSPEC,
				FormFactory.RELATED_GAP_ROWSPEC,
				FormFactory.DEFAULT_ROWSPEC,
				FormFactory.RELATED_GAP_ROWSPEC,
				FormFactory.DEFAULT_ROWSPEC,
				FormFactory.RELATED_GAP_ROWSPEC,
				FormFactory.RELATED_GAP_ROWSPEC,
				FormFactory.RELATED_GAP_ROWSPEC,
				FormFactory.DEFAULT_ROWSPEC,}));
		
		JRadioButton rdbtnNewRadioButton_2 = new JRadioButton("Use non-evolutionary model");
		panel_4.add(rdbtnNewRadioButton_2, "2, 2");
		
		JRadioButton rdbtnNewRadioButton_1 = new JRadioButton("Infer from alignment using FastTree");
		panel_4.add(rdbtnNewRadioButton_1, "2, 4");
		
		JRadioButton rdbtnNewRadioButton = new JRadioButton("Select from file (Newick format):");
		panel_4.add(rdbtnNewRadioButton, "2, 6");
		
		JButton btnBrowse = new JButton("Browse..");
		panel_4.add(btnBrowse, "2, 8, right, default");
		
		textField_1 = new JTextField();
		panel_4.add(textField_1, "4, 8, fill, default");
		textField_1.setColumns(10);
		
		JButton btnViewInFigtree = new JButton("View selected phylogeny in FigTree");
		panel_4.add(btnViewInFigtree, "4, 12, right, default");
		
		JPanel panel_2 = new JPanel();
		panel_6.add(panel_2);
		panel_2.setBorder(new TitledBorder(UIManager.getBorder("TitledBorder.border"), "Run options", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		panel_2.setLayout(new BorderLayout(0, 0));
		
		JPanel panel_5 = new JPanel();
		panel_2.add(panel_5, BorderLayout.CENTER);
		GridBagLayout gbl_panel_5 = new GridBagLayout();
		gbl_panel_5.columnWidths = new int[]{231, 0, 231, 0};
		gbl_panel_5.rowHeights = new int[]{24, 0, 0};
		gbl_panel_5.columnWeights = new double[]{1.0, 0.0, 0.0, Double.MIN_VALUE};
		gbl_panel_5.rowWeights = new double[]{0.0, 0.0, Double.MIN_VALUE};
		panel_5.setLayout(gbl_panel_5);
		
		JButton btnNewButton_1 = new JButton("Run");
		btnNewButton_1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
			}
		});
		
		JLabel lblNumberOfProcessors = new JLabel("Number of processors");
		GridBagConstraints gbc_lblNumberOfProcessors = new GridBagConstraints();
		gbc_lblNumberOfProcessors.anchor = GridBagConstraints.EAST;
		gbc_lblNumberOfProcessors.insets = new Insets(0, 0, 5, 5);
		gbc_lblNumberOfProcessors.gridx = 0;
		gbc_lblNumberOfProcessors.gridy = 0;
		panel_5.add(lblNumberOfProcessors, gbc_lblNumberOfProcessors);
		
		JSpinner spinner = new JSpinner();
		spinner.setModel(new SpinnerNumberModel(new Integer(1), null, null, new Integer(1)));
		GridBagConstraints gbc_spinner = new GridBagConstraints();
		gbc_spinner.insets = new Insets(0, 0, 5, 5);
		gbc_spinner.anchor = GridBagConstraints.WEST;
		gbc_spinner.fill = GridBagConstraints.VERTICAL;
		gbc_spinner.gridx = 1;
		gbc_spinner.gridy = 0;
		panel_5.add(spinner, gbc_spinner);
		GridBagConstraints gbc_btnNewButton_1 = new GridBagConstraints();
		gbc_btnNewButton_1.anchor = GridBagConstraints.EAST;
		gbc_btnNewButton_1.gridx = 2;
		gbc_btnNewButton_1.gridy = 1;
		panel_5.add(btnNewButton_1, gbc_btnNewButton_1);
	}

}
