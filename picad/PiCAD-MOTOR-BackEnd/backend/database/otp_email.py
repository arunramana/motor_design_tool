import smtplib  # Import smtplib for sending email
from email.mime.text import MIMEText  # Import MIMEText for the email body
from email.mime.multipart import MIMEMultipart  # Import MIMEMultipart for the email structure
import random


def send_alert_mail(subject, content, email):
    """Send an alert email with the given subject and content."""
    
    # Check if subject and content are empty
    if not subject.strip() or not content.strip():
        return 1  # Return 1 if either subject or content is empty

    # Prepare email
    email_message = create_email(subject.strip(), content.strip(), email.strip())
    
    # SMTP server configuration
    smtp_settings = {
        "host": "smtpout.secureserver.net",
        "port": 587,
        "username": "tool-info@mojomotor.in",
        "password": "Pilabz@2024"
    }

    # Send the email
    if send_email(email_message, smtp_settings):
        return 0  # Return 0 if the email was sent successfully
    else:
        return 1  # Return 1 if there was an error

def create_email(subject, content, email):
    """Create an email message with the specified subject and content."""
    
    email_message = MIMEMultipart()  # Create a MIMEMultipart object to structure the email
    email_message['From'] = 'tool-info@mojomotor.in'  # Set the sender's email address
    email_message['To'] = email  # Set the recipient's email address
    email_message['Subject'] = subject  # Set the subject of the email
    email_message.attach(MIMEText(content, 'plain'))  # Attach the body of the email as plain text
    
    return email_message

def send_email(email_message, smtp_settings):
    """Send an email using the specified SMTP settings."""
    
    try:
        # Connect to the SMTP server
        with smtplib.SMTP(smtp_settings["host"], smtp_settings["port"]) as server:
            server.starttls()  # Start TLS for security
            server.login(smtp_settings["username"], smtp_settings["password"])  # Log in to the SMTP server
            # Send the email
            server.sendmail(
                email_message['From'],  # From address
                email_message['To'].split(","),  # To and CC addresses as a list
                email_message.as_string()  # Email content as a string
            )
        return True  # Return True if the email was sent successfully
    except Exception as exc:
        print(f"Error sending email: {exc}")  # Print the error if an exception occurs
        return False  # Return False if there was an error


def send_otp_email(email, otp):
    
    # print(email, 1)
    
    subject = "Reset Password"
    content = f"""A request to reset your password or unlock your account was made for your Email ID, {email}. To continue with this request, please enter the code below on the verification page: \n\n {otp} . \n\n NOTE : This OTP is valid for only 10 minutes and please do not reply to this email. \n\n\n\nSincerely,\nMotormojo Team."""
    result = send_alert_mail(subject, content, email)
    
    if result == 0:
        return "Email sent successfully.", 200
    else:
        return "Error sending email.", 500

def signup_email_otp(email, otp):
    
    # print(email, 1)
    
    subject = "Signup Verification"
    content = f"""A request to verify your Email ID, {email}. To continue with this request, please enter the code below on the verification page: \n\n {otp} . \n\n NOTE : This OTP is valid for only 10 minutes and please do not reply to this email. \n\n\n\nSincerely,\nMotormojo Team."""
    result = send_alert_mail(subject, content, email)
    
    if result == 0:
        return "Email sent successfully.", 200
    else:
        return "Error sending email.", 500